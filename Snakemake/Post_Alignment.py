shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for POST Alignment
# snakemake --snakefile Post_Alignment.py --configfile Yoko.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile Post_Alignment.py --configfile Yoko.json --rulegraph | dot -Tsvg > Post_Alignment.svg
# ################################### IMPORT ##################################### #


import os
import sys
import re
sys.path.append(os.path.abspath("./library/"))
import utility
# ################################### FUNCTIONS ################################## #


def build_design_Dict(metadata_Dict):
	#
	design_Dict = {}
	for sample, sample_Dict in metadata_Dict.items():
		#
		if sample_Dict["Design"] not in design_Dict:
			#
			design_Dict[sample_Dict["Design"]] = {}
			design_Dict[sample_Dict["Design"]]["Case"] = []
			design_Dict[sample_Dict["Design"]]["Control"] = []
			if sample_Dict["Type"] == "CASE":
				#
				design_Dict[sample_Dict["Design"]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				design_Dict[sample_Dict["Design"]]["Control"].append(sample)
			else:
				pass
		elif sample_Dict["Design"] in design_Dict:
			#
			if sample_Dict["Type"] == "CASE":
				#
				design_Dict[sample_Dict["Design"]]["Case"].append(sample)
			elif sample_Dict["Type"] == "CONTROL":
				#
				design_Dict[sample_Dict["Design"]]["Control"].append(sample)
			else:
				pass
	return design_Dict
# ################################### CONFIGURATION ############################## #

# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
TITLE = config_general_Dict["TITLE"]
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
LAYOUT = config_data_Dict["LAYOUT"].lower()
GENOME = config_data_Dict["GENOME"].lower()
READ_LENGTH = config_data_Dict["READ_LENGTH"]
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#CONDA
config_conda_Dict = config["CONDA"]
CONDA_INIT = config_conda_Dict["CONDA_INIT"]
CONDA_PY2 = config_conda_Dict["ATAC_Seq_py2"]
CONDA_PY3 = config_conda_Dict["ATAC_Seq_py3"]
ACTIVATE_CONDA_PY2 = config_conda_Dict["ACTIVATE_PY2"]
ACTIVATE_CONDA_PY3 = config_conda_Dict["ACTIVATE_PY3"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#CLUSTER
config_cluster_Dict = config["CLSUTER_CONFIG"]
PROCESSORS = config_cluster_Dict["PROCESSORS"]
MEMORY = config_cluster_Dict["MEMORY"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#METADATA
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
metadata_Dict = utility.build_metadata_dict(METADATA_FILE)
design_Dict = build_design_Dict(metadata_Dict)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#UTILITIES
config_utilities_Dict = config["UTILITIES"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#ALIGNMENT
config_alignment_Dict = config["ALIGNMENT"]
config_qualimap_Dict = config["QUALIMAP"][GENOME]
config_reference_Dict = config["REFERENCE"][GENOME]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#POST_ALIGNMENT
config_post_alignment_Dict = config["POST_ALIGNMENT"]
if LAYOUT == "paired":
	FILTER_MAPQ = config_post_alignment_Dict["FILTER_MAPQ_PAIRED"]
else:
	FILTER_MAPQ = config_post_alignment_Dict["FILTER_MAPQ_SINGLE"]
CHROMOSOME_FILTER_PROCESS = utility.build_snakemake_awk(config_post_alignment_Dict["FILTER_CHROMOSOME"])
# ------------------------------------
# ################################### WILDCARDS ################################ #

alignment_List = []
post_alignment_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		post_alignment_List
# ################################### PIPELINE RULES ########################## #

#+++++++++++++++++++++++++++++
##REALM3: POST-ALIGNMENT
#+++++++++++++++++++++++++++++

rule Post_Alignment:
	input:
		bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
		bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
	output:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Post_Alignment: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			#
			module load samtools/1.9
			module load bedtools/2.27.1
			module load deeptools/3.1.3
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/post_alignment
			mkdir -p $QC_PATH
			#
			#
			if [ ! -f ./Script/{GENOME}.blacklist.bed.gz ]; then
				wget {config_reference_Dict[BLACK_LIST]} -O ./Script/{GENOME}.blacklist.bed.gz
			fi
			AWK_COMMAND='{config_post_alignment_Dict[FILTER_CHROMOSOME]}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "deeptools/3.1.3" | tee >(cat >&2)
			printf "%s\\n" "Filtering bam file, correctGC Bias, drop short alignment"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.bam_index}"  | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.processed_bam}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.processed_bam_index}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "samtools view --threads {threads} -h {FILTER_MAPQ} {input.bam} | awk -F'\\t' '$AWK_COMMAND' | samtools view --threads {threads} -Shb - > {output.processed_bam}.filt" | tee >(cat >&2)
			printf "%s\\n" "bedtools intersect -v -abam {output.processed_bam}.filt -b <(zcat -f ./Script/{GENOME}.blacklist.bed.gz ) > {output.processed_bam}.blk.filt" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.blk.filt -o {output.processed_bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} -b {output.processed_bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "alignmentSieve --bam {output.processed_bam}.tmp --numberOfProcessors {threads} --minFragmentLength {READ_LENGTH} --outFile {output.processed_bam}.fraglen --filterMetrics $QC_PATH/{wildcards.sample}.alignmentSieve.txt"  | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.fraglen -o {output.processed_bam}" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} -b {output.processed_bam}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.processed_bam}.filt {output.processed_bam}.blk.filt {output.processed_bam}.tmp {output.processed_bam}.tmp.bai {output.processed_bam}.fraglen" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			samtools view --threads {threads} -h {FILTER_MAPQ} {input.bam} | awk -F'\\t' '$AWK_COMMAND' | samtools view --threads {threads} -Shb - > {output.processed_bam}.filt
			bedtools intersect -v -abam {output.processed_bam}.filt -b <(zcat -f ./Script/{GENOME}.blacklist.bed.gz ) > {output.processed_bam}.blk.filt
			samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.blk.filt -o {output.processed_bam}.tmp
			samtools index -@ {threads} -b {output.processed_bam}.tmp
			alignmentSieve --bam {output.processed_bam}.tmp --numberOfProcessors {threads} --minFragmentLength {READ_LENGTH} --outFile {output.processed_bam}.fraglen --filterMetrics $QC_PATH/{wildcards.sample}.alignmentSieve.txt
			samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.fraglen -o {output.processed_bam}
			samtools index -@ {threads} -b {output.processed_bam}
			rm -rf {output.processed_bam}.filt {output.processed_bam}.blk.filt {output.processed_bam}.tmp {output.processed_bam}.tmp.bai {output.processed_bam}.fraglen
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
		""")
