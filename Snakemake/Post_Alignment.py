shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for POST Alignment
# snakemake --snakefile Post_Alignment.py --configfile Yoko.json --cores=50 -j 10 --local-cores=10
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
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{sample}.samtools.txt".format(design=sample_Dict["Design"], sample=sample))
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{sample}.processed.samtools.txt".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{pooled_case}.samtools.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{pooled_case}.processed.samtools.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
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
		bowtie2_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
		bowtie2_index_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
	output:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_index_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Post_Alignment: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			########################################################
			#
			if [ ! -f ./Script/{GENOME}.blacklist.bed.gz ]; then
				wget {config_reference_Dict[BLACK_LIST]} -O ./Script/{GENOME}.blacklist.bed.gz
			fi
			#
			AWK_COMMAND='{config_post_alignment_Dict[FILTER_CHROMOSOME]}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "GENERAL FILTERING" | tee >(cat >&2)
			printf "%s\\n" "samtools view --threads {threads} -h {FILTER_MAPQ} {input.bowtie2_bam} | awk -F'\\t' '$AWK_COMMAND' | samtools view --threads {threads} -Shb - > {output.processed_bam}.filt" | tee >(cat >&2)
			printf "%s\\n" "bedtools intersect -v -abam {output.processed_bam}.filt -b <(zcat -f ./Script/{GENOME}.blacklist.bed.gz ) > {output.processed_bam}.blk.filt" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.blk.filt -o {output.processed_bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} -b {output.processed_bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			samtools view --threads {threads} -h {FILTER_MAPQ} {input.bowtie2_bam} | awk -F'\\t' '$AWK_COMMAND' | samtools view --threads {threads} -Shb - > {output.processed_bam}.filt
			bedtools intersect -v -abam {output.processed_bam}.filt -b <(zcat -f ./Script/{GENOME}.blacklist.bed.gz ) > {output.processed_bam}.blk.filt
			samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.blk.filt -o {output.processed_bam}.tmp
			samtools index -@ {threads} -b {output.processed_bam}.tmp
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "FRAGMENT LENGTH FILTER" | tee >(cat >&2)
			printf "%s\\n" "alignmentSieve --bam {output.processed_bam}.tmp --numberOfProcessors {threads} --minFragmentLength 75 --outFile {output.processed_bam}.fraglen"  | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.fraglen -o {output.processed_bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} -b {output.processed_bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			alignmentSieve --bam {input.bowtie2_bam} --numberOfProcessors {threads} --minFragmentLength 75 --outFile {output.processed_bam}.fraglen --filteredOutReads {output.processed_bam}.discarded
			samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.fraglen -o {output.processed_bam}.tmp
			samtools index -@ {threads} -b {output.processed_bam}.tmp
			samtools sort --threads {threads} -m 2G -O bam {output.processed_bam}.discarded -o {output.processed_bam}.discarded.bam
			samtools index -@ {threads} -b {output.processed_bam}.discarded.bam
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "CORRECT GC BIAS" | tee >(cat >&2)
			printf "%s\\n" "computeGCBias --bamfile {output.processed_bam}.tmp --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} --genome {config_reference_Dict[2BIT]} --fragmentLength 150 --GCbiasFrequenciesFile {output.processed_bam}.GC_freq" | tee >(cat >&2)
			printf "%s\\n" "correctGCBias --bamfile {output.processed_bam}.tmp --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} --genome {config_reference_Dict[2BIT]} --GCbiasFrequenciesFile {output.processed_bam}.GC_freq --correctedFile {output.processed_bam}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			computeGCBias --bamfile {output.processed_bam}.tmp --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} --genome {config_reference_Dict[2BIT]} --fragmentLength 150 --GCbiasFrequenciesFile {output.processed_bam}.GC_freq
			correctGCBias --bamfile {output.processed_bam}.tmp --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} --genome {config_reference_Dict[2BIT]} --GCbiasFrequenciesFile {output.processed_bam}.GC_freq --correctedFile {output.processed_bam}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			
		""")


rule Post_Alignment_QC:
	input:
		processed_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
	output:
		flagstats = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{sample}.processed.samtools.txt",
		spp_phantomPeak = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{sample}.processed.spp.pdf",
		spp_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{sample}.processed.spp.txt",
		complexity_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{sample}.processed.complexity.txt",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Post_Alignment_QC: {wildcards.design}|{wildcards.sample}"
	run:
		
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/post_alignment
			mkdir -p $QC_PATH
			unset DISPLAY
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "qualimap bamqc -bam {input.processed_Bam} -nt {threads} -outdir $QC_PATH/{wildcards.sample}_processed_qualimap {config_qualimap_Dict}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			qualimap bamqc -bam {input.processed_Bam} -nt {threads} -outdir $QC_PATH/{wildcards.sample}_processed_qualimap {config_qualimap_Dict}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "fastqc -o $QC_PATH --format bam --threads {threads} {input.processed_Bam}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			fastqc -o $QC_PATH --format bam --threads {threads} {input.processed_Bam}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools flagstat --threads {threads} {input.processed_Bam} > $QC_PATH/{wildcards.sample}.processed.samtools.txt" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			samtools flagstat --threads {threads} {input.processed_Bam} > $QC_PATH/{wildcards.sample}.processed.samtools.txt
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			#
			if [ ! -f ./R_Script/run_spp_nodups.R ]; then
				wget {config_utilities_Dict[SPP]} -O ./R_Script/run_spp_nodups.R
			fi
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "Rscript ./R_Script/run_spp_nodups.R -c={input.processed_Bam} -p={threads} -savp={output.spp_phantomPeak} -speak=0 -rf > {output.spp_report}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			Rscript ./R_Script/run_spp_nodups.R -c={input.processed_Bam} -p={threads} -savp={output.spp_phantomPeak} -speak=0 -rf > {output.spp_report}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "python ./Python_Script/bamChipQC.py -i {input.processed_Bam} -c {threads} -o {output.complexity_report}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			python ./Python_Script/bamChipQC.py -i {input.processed_Bam} -c {threads} -o {output.complexity_report}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			
		""")