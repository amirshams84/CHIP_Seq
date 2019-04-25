# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for alignment
# snakemake --snakefile alignment.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile alignment.py --configfile Encode.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
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


def get_case_bam(wildcards):
	"""
	"""
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{case}.bam".format(design=wildcards.design, case=sample))
	return bam_List
# ################################### CONFIGURATION ############################## #

# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
TITLE = config_general_Dict["TITLE"]
EXECUTION_MODE = config_general_Dict["EXECUTION_MODE"]
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
# ################################### WILDCARDS ################################ #


pre_process_List = []
alignment_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	pre_process_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq".format(design=sample_Dict["Design"], sample=sample))
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		alignment_List
# ################################### PIPELINE RULES ########################## #


if LAYOUT == "paired":
	rule Alignment_Paired:
		"""
		"""
		input:
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R2.processed.fastq",
		output:
			bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
			bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "Alignment_Paired: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fwd_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq", "", each_fastq_basename)
			shell("""
				#
				module load bowtie/2-2.3.5
				module load samtools/1.9
				module load picard/2.18.27
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				#
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "bowtie/2-2.3.5" | tee >(cat >&2)
				printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
				printf "%s\\n" "Align processed reads to specified genome and sort it, mark it with duplicate, index it"  | tee >(cat >&2)
				printf "INPUT1: %s\\n" "{input.processed_fwd_fastq}"  | tee >(cat >&2)
				printf "INPUT2: %s\\n" "{input.processed_rev_fastq}"  | tee >(cat >&2)
				printf "OUTPUT1: %s\\n" "{output.bam}"  | tee >(cat >&2)
				printf "OUTPUT2: %s\\n" "{output.bam_index}"  | tee >(cat >&2)
				printf "OUTPUT3: %s\\n" "$QC_PATH/{each_fastq_begining}.alignment.txt"  | tee >(cat >&2)
				printf "OUTPUT3: %s\\n" "$QC_PATH/{each_fastq_begining}.picard.txt"  | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} 2> $QC_PATH/{each_fastq_begining}.alignment.txt | samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bam}.tmp - " | tee >(cat >&2)
				printf "%s\\n" "java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp OUTPUT={output.bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}" | tee >(cat >&2)
				printf "%s\\n" "samtools index -@ {threads} -b {output.bam}" | tee >(cat >&2)
				printf "%s\\n" "rm -rf {output.bam}.tmp" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} 2> $QC_PATH/{each_fastq_begining}.alignment.txt | samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bam}.tmp -
				java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp OUTPUT={output.bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}
				samtools index -@ {threads} -b {output.bam}
				rm -rf {output.bam}.tmp
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			""")

elif LAYOUT == "single":
	rule Alignment_Single:
		"""
		"""
		input:
			processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq",
		output:
			bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
			bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "Alignment_Single: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq", "", each_fastq_basename)
			shell("""
				#
				module load bowtie/2-2.3.5
				module load samtools/1.9
				module load picard/2.18.27
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				#
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "bowtie/2-2.3.5" | tee >(cat >&2)
				printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
				printf "%s\\n" "Align processed reads to specified genome and sort it, mark it with duplicate, index it"  | tee >(cat >&2)
				printf "INPUT1: %s\\n" "{input.processed_fwd_fastq}"  | tee >(cat >&2)
				printf "INPUT2: %s\\n" "{input.processed_rev_fastq}"  | tee >(cat >&2)
				printf "OUTPUT1: %s\\n" "{output.bam}"  | tee >(cat >&2)
				printf "OUTPUT2: %s\\n" "{output.bam_index}"  | tee >(cat >&2)
				printf "OUTPUT3: %s\\n" "$QC_PATH/{each_fastq_begining}.alignment.txt"  | tee >(cat >&2)
				printf "OUTPUT4: %s\\n" "$QC_PATH/{each_fastq_begining}.picard.txt"  | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -U {input.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.alignment.txt | samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bam}.tmp - " | tee >(cat >&2)
				printf "%s\\n" "java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp OUTPUT={output.bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}" | tee >(cat >&2)
				printf "%s\\n" "samtools index -@ {threads} -b {output.bam}" | tee >(cat >&2)
				printf "%s\\n" "rm -rf {output.bam}.tmp" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} 2> $QC_PATH/{each_fastq_begining}.alignment.txt | samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bam}.tmp -
				java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.bam}.tmp OUTPUT={output.bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}
				samtools index -@ {threads} -b {output.bam}
				rm -rf {output.bam}.tmp
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			""")


rule Pooling_Case_Replicates:
	input:
		bam_List = get_case_bam
	output:
		pooled_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam",
		pooled_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam.bai",
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Pooling_Case_Replicates: {wildcards.design}|{wildcards.pooled_case}"
	run:
		shell("""
			#
			module load samtools/1.9
			module load picard/2.18.27
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
			mkdir -p $QC_PATH
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "Merging Biological replicates bam files into single pooled sample, sort mark duplicates, and index"  | tee >(cat >&2)
			declare -a bam_List=({input.bam_List})
			for case_index in "${{!bam_List[@]}}"
			do
				index=$(($case_index+1))
				printf "INPUT%s: %s\\n" "${{index}}:" "${{bam_List[$case_index]}}" | tee >(cat >&2)
			done
			printf "OUTPUT1: %s\\n" "{output.pooled_bam}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.pooled_bam_index}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "samtools merge --threads {threads} {output.pooled_bam} {input.bam_List}" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -m 2G -O bam {output.pooled_bam}.unsorted -o {output.pooled_bam}.sorted" | tee >(cat >&2)
			printf "%s\\n" "java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.pooled_bam}.sorted OUTPUT={output.pooled_bam} METRICS_FILE=$QC_PATH/{wildcards.pooled_case}.picard.txt {config_alignment_Dict[PICARD]}" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} {output.pooled_bam}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.pooled_bam}.unsorted" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			samtools merge --threads {threads} {output.pooled_bam}.unsorted {input.bam_List}
			samtools sort --threads {threads} -m 2G -O bam {output.pooled_bam}.unsorted -o {output.pooled_bam}.sorted
			java -Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={output.pooled_bam}.sorted OUTPUT={output.pooled_bam} METRICS_FILE=$QC_PATH/{wildcards.pooled_case}.picard.txt {config_alignment_Dict[PICARD]}
			samtools index -@ {threads} {output.pooled_bam}
			rm -rf {output.pooled_bam}.unsorted {output.pooled_bam}.sorted
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
		""")
