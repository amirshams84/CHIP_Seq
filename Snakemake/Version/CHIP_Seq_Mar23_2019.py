shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-19-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for CHIP_Seq data analysis
# snakemake --snakefile CHIP_Seq.py --configfile Encode.json --cores=50 -j 1 --local-cores=10
# snakemake --snakefile CHIP_Seq.py --configfile Yoko.json --cores=50 -j 1 --local-cores=10
# snakemake --snakefile CHIP_Seq.py --configfile CHIP_Seq.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
# snakemake --snakefile CHIP_Seq.py --configfile Yoko.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
# ################################### IMPORT ##################################### #
#
import os
import re
import sys
import pickle
import glob
import collections
import pandas
from shutil import copyfile
from importlib import reload
sys.path.append(os.path.abspath("./Python_Script/"))
import utility
# ################################### FUNCTIONS ################################## #


def get_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + SAMPLE_DELIMITER + "*" + "." + DATA_SUFFIX, recursive=True)


def get_forward_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + SAMPLE_DELIMITER + "*" + FORWARD_DELIMITER + "*" + "." + DATA_SUFFIX, recursive=True)


def get_case_bam(wildcards):
	"""
	"""
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			if sample_Dict["Type"] == "CASE":
				bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{case}.bam".format(design=wildcards.design, case=sample))
	return bam_List


def get_bam(wildcards):
	"""
	"""
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=wildcards.design, sample=sample))
	return bam_List


def get_pooled_bam(wildcards):
	"""
	"""
	pooled_bam_List = []
	pooled_bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
	return pooled_bam_List


def get_bigwig(wildcards):
	"""
	"""
	bigwig_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			bigwig_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample}.bigwig".format(design=wildcards.design, sample=sample))
	return bigwig_List


def get_pooled_bigwig(wildcards):
	"""
	"""
	pooled_bigwig_List = []
	pooled_bigwig_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{pooled_case}.bigwig".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
	return pooled_bigwig_List


def get_processed_bam(wildcards):
	"""
	"""
	bam_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=wildcards.design, sample=sample))
	return bam_List


def get_pooled_processed_bam(wildcards):
	"""
	"""
	pooled_bam_List = []
	pooled_bam_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
	return pooled_bam_List


def get_processed_bigwig(wildcards):
	"""
	"""
	bigwig_List = []
	for sample, sample_Dict in metadata_Dict.items():
		if sample_Dict["Design"] == wildcards.design:
			bigwig_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bigwig".format(design=wildcards.design, sample=sample))
	return bigwig_List


def get_pooled_processed_bigwig(wildcards):
	"""
	"""
	pooled_bigwig_List = []
	pooled_bigwig_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bigwig".format(design=wildcards.design, pooled_case="_POOLED_".join(design_Dict[wildcards.design]["Case"])))
	return pooled_bigwig_List

# ################################### CONFIGURATION ############################## #

# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
INFOLINK = config_general_Dict["INFOLINK"]
TITLE = config_general_Dict["TITLE"]
GENOME = config_general_Dict["GENOME"]
GENOME = GENOME.lower()
SAMPLE_DELIMITER = config_general_Dict["SAMPLE_DELIMITER"]
FORWARD_DELIMITER = config_general_Dict["FORWARD_DELIMITER"]
REVERSE_DELIMITER = config_general_Dict["REVERSE_DELIMITER"]
FORMAT = config_general_Dict["FORMAT"]
LAYOUT = config_general_Dict["LAYOUT"].lower()
PLATFORM = config_general_Dict["PLATFORM"]
DATA_SUFFIX = config_general_Dict["DATA_SUFFIX"].lower()
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
DATADIR = utility.fix_path(config_general_Dict["DATADIR"])
# ++++++++++++++++++++++++++++++++++++
#
if FORMAT.lower() != "fastq":
	print("This pipeline only can use fastq format")
	print("Aborting!!!!")
	sys.exit(2)
else:
	pass

if DATA_SUFFIX[0] == ".":
	DATA_SUFFIX = DATA_SUFFIX[1:]
if DATA_SUFFIX.lower() not in ["fastq", "fastq.gz"]:
	print("This pipeline only can use fastq or fastq.gz data")
	print("Aborting!!!!")
	sys.exit(2)
else:
	pass
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
#UTILITIES
config_utilities_Dict = config["UTILITIES"]
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
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#PRE_PROCESS
config_pre_process_Dict = config["PRE_PROCESS"]
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
# ++++++++++++++++++++++++++++++++++++
#PEAK_CALLING
config_peak_calling_Dict = config["PEAK_CALLING"][GENOME]
if LAYOUT == "paired":
	MACS2_NARROW_PARAMETERS = config_peak_calling_Dict["MACS2_NARROW_PAIRED"]
	MACS2_BROAD_PARAMETERS = config_peak_calling_Dict["MACS2_BROAD_PAIRED"]
else:
	MACS2_NARROW_PARAMETERS = config_peak_calling_Dict["MACS2_NARROW_SINGLE"]
	MACS2_BROAD_PARAMETERS = config_peak_calling_Dict["MACS2_BROAD_SINGLE"]
# ------------------------------------
# ################################### WILDCARDS ################################ #
design_Dict = {}
pre_process_List = []
alignment_List = []
post_alignment_List = []
deeptools_List = []
peak_calling_List = []

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
	#
	#PRE_PROCESS
	##PROCESSED_FASTQ
	pre_process_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	#ALIGNMENT
	##BOWTIE2
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	##BAM_QC
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/alignment/{sample}.samtools.txt".format(design=sample_Dict["Design"], sample=sample))
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	##POST_ALIGNMENT_QC
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/post_alignment/{sample}.processed.samtools.txt".format(design=sample_Dict["Design"], sample=sample))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))

###############

cross_correlate_List = []

for design in design_Dict:
	#
	##POOLING
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/alignment/{pooled_case}.samtools.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/post_alignment/{pooled_case}.processed.samtools.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	
	#post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/correlation/{design}.read_coverage.txt".format(design=design))
	##VERSUS
	for case in design_Dict[design]["Case"]:
		for control in design_Dict[design]["Control"]:
			#
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
	##
	for control in design_Dict[design]["Control"]:
		#
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{pooled_case}_VS_{control}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
	#
	##BAM_CORRELATE
	cross_correlate_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/post_alignment/{design}.read_coverage.txt".format(design=design))
# ################################### PIPELINE FLOW ############################## #


rule End_Point:
	input:
		pre_process_List
		+ alignment_List
		+ post_alignment_List
		+ peak_calling_List
# ################################### PIPELINE RULES ############################ #


#+++++++++++++++++++++++++++++
##REALM1: PRE_PROCESS
#+++++++++++++++++++++++++++++
if LAYOUT == "paired":
	#
	rule PreProcess_Paired:
		input:
			fastq_List = get_forward_fastq
		output:
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz",
		threads: PROCESSORS
		message: "PreProcess: {wildcards.design}|{wildcards.sample}"
		resources:
			mem_mb = MEMORY
		run:
			for each_fastq in input.fastq_List:
				#
				each_fastq_reverse = each_fastq.replace(FORWARD_DELIMITER, REVERSE_DELIMITER)
				each_fastq_basename = os.path.basename(each_fastq)
				each_fastq_begining = re.sub("." + DATA_SUFFIX, "", each_fastq_basename)
				shell("""
					{ACTIVATE_CONDA_PY3}
					QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/report/pre_process
					mkdir -p $QC_PATH
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "cutadapt {config_pre_process_Dict[CUTADAPT_PAIRED]} --cores={threads} --output={output.processed_fwd_fastq}.tmp --paired-output={output.processed_rev_fastq}.tmp {each_fastq} {each_fastq_reverse} > $QC_PATH/{each_fastq_begining}.cutadapt.txt" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					start_time="$(date -u +%s)"
					#
					##
					cutadapt {config_pre_process_Dict[CUTADAPT_PAIRED]} --cores={threads} --output={output.processed_fwd_fastq}.tmp --paired-output={output.processed_rev_fastq}.tmp {each_fastq} {each_fastq_reverse} > $QC_PATH/{each_fastq_begining}.cutadapt.txt
					##
					#
					end_time="$(date -u +%s)"
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					cat {output.processed_fwd_fastq}.tmp >> {output.processed_fwd_fastq}
					cat {output.processed_rev_fastq}.tmp >> {output.processed_rev_fastq}
					rm {output.processed_fwd_fastq}.tmp {output.processed_rev_fastq}.tmp
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "fastqc -o $QC_PATH --format fastq --threads {threads} {each_fastq}" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					start_time="$(date -u +%s)"
					#
					##
					fastqc -o $QC_PATH --format fastq --threads {threads} {each_fastq}
					##
					#
					end_time="$(date -u +%s)"
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "fastqc -o $QC_PATH --format fastq --threads {threads} {each_fastq_reverse}" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					start_time="$(date -u +%s)"
					#
					##
					fastqc -o $QC_PATH --format fastq --threads {threads} {each_fastq_reverse}
					##
					#
					end_time="$(date -u +%s)"
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				""")
else:
	#
	rule PreProcess:
		input:
			fastq_List = get_fastq
		output:
			processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		threads: PROCESSORS
		message: "PreProcess: {wildcards.design}|{wildcards.sample}"
		resources:
			mem_mb = MEMORY
		run:
			for each_fastq in input.fastq_List:
				#
				each_fastq_basename = os.path.basename(each_fastq)
				each_fastq_begining = re.sub("." + DATA_SUFFIX, "", each_fastq_basename)
				shell("""
					{ACTIVATE_CONDA_PY3}
					QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/report/pre_process
					mkdir -p $QC_PATH
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "cutadapt {config_pre_process_Dict[CUTADAPT_SINGLE]} --cores={threads} {each_fastq} >> {output.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.cutadapt.txt"  | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					start_time="$(date -u +%s)"
					#
					##
					cutadapt {config_pre_process_Dict[CUTADAPT_SINGLE]} --cores={threads} {each_fastq} >> {output.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.cutadapt.txt
					##
					#
					end_time="$(date -u +%s)"
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "fastqc -o $QC_PATH --format fastq --threads {threads} {each_fastq}" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					start_time="$(date -u +%s)"
					#
					##
					fastqc -o $QC_PATH --format fastq --threads {threads} {each_fastq}
					##
					#
					end_time="$(date -u +%s)"
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				""")


#+++++++++++++++++++++++++++++
##REALM2: ALIGNMENT
#+++++++++++++++++++++++++++++
if LAYOUT == "paired":
	#
	rule Alignment_Paired:
		input:
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R2.processed.fastq.gz",
		output:
			bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
			bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam.bai",
			bowtie2_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.bigwig"
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "Alignment_Paired: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fwd_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq.gz", "", each_fastq_basename)
			shell("""
				{ACTIVATE_CONDA_PY2}
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} 2> $QC_PATH/{each_fastq_begining}.bowtie2.txt | \
				samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bowtie2_Bam}.tmp - " | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				bowtie2 {config_alignment_Dict[BOWTIE2_PAIRED]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -1 {input.processed_fwd_fastq} -2 {input.processed_rev_fastq} 2> $QC_PATH/{each_fastq_begining}.bowtie2.txt | \
				samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bowtie2_Bam}.tmp -
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "export _JAVA_OPTIONS='-Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads}'" | tee >(cat >&2)
				printf "%s\\n" "picard MarkDuplicates INPUT={output.bowtie2_Bam}.tmp OUTPUT={output.bowtie2_Bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				export _JAVA_OPTIONS='-Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads}'
				picard MarkDuplicates INPUT={output.bowtie2_Bam}.tmp OUTPUT={output.bowtie2_Bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "samtools index -@ {threads} -b {output.bowtie2_Bam}" | tee >(cat >&2)
				printf "%s\\n" "rm -rf {output.bowtie2_Bam}.tmp" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				samtools index -@ {threads} -b {output.bowtie2_Bam}
				rm -rf {output.bowtie2_Bam}.tmp
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "bamCoverage --bam {output.bowtie2_Bam} --outFileName {output.bowtie2_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				bamCoverage --bam {output.bowtie2_Bam} --outFileName {output.bowtie2_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			""")
else:
	#
	rule Alignment:
		input:
			processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/pre_process/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.R1.processed.fastq.gz",
		output:
			bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
			bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bam.bai",
			bowtie2_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_POOLED_.*|.*_OVERLAPPED_.*|\\.).)*}.bigwig",
			
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "Alignment: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq.gz", "", each_fastq_basename)
			shell("""
				{ACTIVATE_CONDA_PY2}
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2_SINGLE]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -U {input.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.bowtie2.txt | \
				samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bowtie2_Bam}.tmp - " | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				bowtie2 {config_alignment_Dict[BOWTIE2_SINGLE]} --threads {threads} -x {config_reference_Dict[BOWTIE2_INDEX]} -U {input.processed_fastq} 2> $QC_PATH/{each_fastq_begining}.bowtie2.txt | \
				samtools sort --threads {threads} -m 2G -O bam -T {WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/alignment/{wildcards.sample} -o {output.bowtie2_Bam}.tmp -
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "export _JAVA_OPTIONS='-Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads}'" | tee >(cat >&2)
				printf "%s\\n" "picard MarkDuplicates INPUT={output.bowtie2_Bam}.tmp OUTPUT={output.bowtie2_Bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				export _JAVA_OPTIONS='-Xms10000M -Xmx10000M -XX:ParallelGCThreads={threads}'
				picard MarkDuplicates INPUT={output.bowtie2_Bam}.tmp OUTPUT={output.bowtie2_Bam} METRICS_FILE=$QC_PATH/{each_fastq_begining}.picard.txt {config_alignment_Dict[PICARD]}
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "samtools index -@ {threads} -b {output.bowtie2_Bam}" | tee >(cat >&2)
				printf "%s\\n" "rm -rf {output.bowtie2_Bam}.tmp" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				samtools index -@ {threads} -b {output.bowtie2_Bam}
				rm -rf {output.bowtie2_Bam}.tmp
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "bamCoverage --bam {output.bowtie2_Bam} --outFileName {output.bowtie2_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				bamCoverage --bam {output.bowtie2_Bam} --outFileName {output.bowtie2_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			""")


rule Pooling_Case_Replicates:
	input:
		bam_List = get_case_bam
	output:
		pooled_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam",
		pooled_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam.bai",
		pooled_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bigwig"
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Pooling_Case_Replicates: {wildcards.design}|{wildcards.pooled_case}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools merge --threads {threads} {output.pooled_Bam} {input.bam_List}" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -m 2G -O bam {output.pooled_Bam}.unsorted -o {output.pooled_Bam}" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} {output.pooled_Bam}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.pooled_Bam}.unsorted" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			samtools merge --threads {threads} {output.pooled_Bam}.unsorted {input.bam_List}
			samtools sort --threads {threads} -m 2G -O bam {output.pooled_Bam}.unsorted -o {output.pooled_Bam}
			samtools index -@ {threads} {output.pooled_Bam}
			rm -rf {output.pooled_Bam}.unsorted
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "bamCoverage --bam {output.pooled_Bam} --outFileName {output.pooled_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			bamCoverage --bam {output.pooled_Bam} --outFileName {output.pooled_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
		""")


rule Alignment_QC:
	input:
		bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
		bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.bam.bai",
	output:
		flagstats = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.samtools.txt",
		spp_phantomPeak = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.spp.pdf",
		spp_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.spp.txt",
		complexity_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.complexity.txt",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "ALIGNMENT_QC: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/report/alignment
			mkdir -p $QC_PATH
			unset DISPLAY
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "qualimap bamqc -bam {input.bowtie2_Bam} -nt {threads} -outdir $QC_PATH/{wildcards.sample}_qualimap {config_qualimap_Dict}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			qualimap bamqc -bam {input.bowtie2_Bam} -nt {threads} -outdir $QC_PATH/{wildcards.sample}_qualimap {config_qualimap_Dict}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "fastqc -o $QC_PATH --format bam --threads {threads} {input.bowtie2_Bam}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			fastqc -o $QC_PATH --format bam --threads {threads} {input.bowtie2_Bam}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools flagstat --threads {threads} {input.bowtie2_Bam} > $QC_PATH/{wildcards.sample}.samtools.txt" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			samtools flagstat --threads {threads} {input.bowtie2_Bam} > $QC_PATH/{wildcards.sample}.samtools.txt
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
			printf "%s\\n" "Rscript ./R_Script/run_spp_nodups.R -c={input.bowtie2_Bam} -p={threads} -savp={output.spp_phantomPeak} -speak=0 -rf > {output.spp_report}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			Rscript ./R_Script/run_spp_nodups.R -c={input.bowtie2_Bam} -p={threads} -savp={output.spp_phantomPeak} -speak=0 -rf > {output.spp_report}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "python ./Python_Script/bamChipQC.py -i {input.bowtie2_Bam} -c {threads} -o {output.complexity_report}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			python ./Python_Script/bamChipQC.py -i {input.bowtie2_Bam} -c {threads} -o {output.complexity_report}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
		""")


#+++++++++++++++++++++++++++++
##REALM3: POST-ALIGNMENT
#+++++++++++++++++++++++++++++


rule Post_Alignment:
	input:
		bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.bam",
		bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.bam.bai",
	output:
		processed_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam",
		processed_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam.bai",
		processed_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bigwig"
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Post_Alignment: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}

			#
			if [ ! -f ./Script/{GENOME}.blacklist.bed.gz ]; then
				wget {config_reference_Dict[BLACK_LIST]} -O ./Script/{GENOME}.blacklist.bed.gz
			fi
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			AWK_COMMAND='{config_post_alignment_Dict[FILTER_CHROMOSOME]}'
			printf "%s\\n" "samtools view --threads {threads} -h {FILTER_MAPQ} {input.bowtie2_Bam} | awk -F'\\t' '$AWK_COMMAND' | samtools view --threads {threads} -Shb - > {output.processed_Bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedtools intersect -v -abam {output.processed_Bam}.tmp -b <(zcat -f ./Script/{GENOME}.blacklist.bed.gz ) > {output.processed_Bam}.unsorted" | tee >(cat >&2)
			printf "%s\\n" "samtools sort --threads {threads} -m 2G -O bam {output.processed_Bam}.unsorted -o {output.processed_Bam}" | tee >(cat >&2)
			printf "%s\\n" "samtools index -@ {threads} -b {output.processed_Bam}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.processed_Bam}.unsorted {output.processed_Bam}.tmp" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			samtools view --threads {threads} -h {FILTER_MAPQ} {input.bowtie2_Bam} | awk -F'\\t' '{CHROMOSOME_FILTER_PROCESS}' | samtools view --threads {threads} -Shb - > {output.processed_Bam}.tmp
			bedtools intersect -v -abam {output.processed_Bam}.tmp -b <(zcat -f ./Script/{GENOME}.blacklist.bed.gz ) > {output.processed_Bam}.unsorted
			samtools sort --threads {threads} -m 2G -O bam {output.processed_Bam}.unsorted -o {output.processed_Bam}
			samtools index -@ {threads} -b {output.processed_Bam}
			rm -rf {output.processed_Bam}.unsorted {output.processed_Bam}.tmp
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "bamCoverage --bam {output.processed_Bam} --outFileName {output.processed_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			bamCoverage --bam {output.processed_Bam} --outFileName {output.processed_bigwig} --numberOfProcessors {threads} --effectiveGenomeSize {config_reference_Dict[EFFECTIVE_GENOME_SIZE]} {config_alignment_Dict[DEEPTOOLS]}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			
		""")


rule Post_Alignment_QC:
	input:
		processed_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam",
		processed_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam.bai",
	output:
		flagstats = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.samtools.txt",
		spp_phantomPeak = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.spp.pdf",
		spp_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.spp.txt",
		complexity_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/report/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.complexity.txt",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Post_Alignment_QC: {wildcards.design}|{wildcards.sample}"
	run:
		
		shell("""
			{ACTIVATE_CONDA_PY2}
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/report/post_alignment
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

#+++++++++++++++++++++++++++++
##REALM4: PEAK-CALLING
#+++++++++++++++++++++++++++++

rule Peak_Calling_Narrow:
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam",
		processed_index_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam.bai",
	output:
		#
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.gz",
		narrowPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.bb",
		narrowPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.FE.bigwig",
		narrowPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.PV.bigwig",
		narrowPeak_QV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.narrowPeak.QV.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Narrow: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {input.processed_bam})
			sample_Name=${{sample_Name%.processed.bam}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling
			#
			if [ ! -f ./Script/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O ./Script/bigNarrowPeak.as
			fi
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bam} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_bam} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			printf "%s\\n" "awk '$AWK_COMMAND' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.narrowPeak_bed}.sorted {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted
			awk 'BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			rm -rf {output.narrowPeak_bed}.sorted {output.narrowPeak_bed}.tmp
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_FE_bigwig}.bdg.tmp > {output.narrowPeak_FE_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_FE_bigwig}.bdg.tmp > {output.narrowPeak_FE_bigwig}.bdg
			bedGraphToBigWig {output.narrowPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method ppois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method ppois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_PV_bigwig}.bdg.tmp > {output.narrowPeak_PV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_PV_bigwig}.bdg.tmp > {output.narrowPeak_PV_bigwig}.bdg
			bedGraphToBigWig {output.narrowPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method qpois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method qpois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_qpois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_QV_bigwig}.bdg.tmp > {output.narrowPeak_QV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_QV_bigwig}.bdg.tmp > {output.narrowPeak_QV_bigwig}.bdg
			bedGraphToBigWig {output.narrowPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			if [ -f {output.narrowPeak_QV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_qpois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_qpois.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			fi
		""")


rule Peak_Calling_Narrow_Controlled:
	input:
		processed_case_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
	output:
		#
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.gz",
		narrowPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.bb",
		narrowPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.FE.bigwig",
		narrowPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.PV.bigwig",
		narrowPeak_QV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.QV.bigwig",
		
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Narrow_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {output.narrowPeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling
			#
			if [ ! -f ./Script/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O ./Script/bigNarrowPeak.as
			fi
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_Bam} --control {input.processed_control_Bam} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_case_Bam} --control {input.processed_control_Bam} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			printf "%s\\n" "awk '$AWK_COMMAND' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.narrowPeak_bed}.sorted {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted
			awk 'BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			rm -rf {output.narrowPeak_bed}.sorted {output.narrowPeak_bed}.tmp
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_FE_bigwig}.bdg.tmp > {output.narrowPeak_FE_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_FE_bigwig}.bdg.tmp > {output.narrowPeak_FE_bigwig}.bdg
			bedGraphToBigWig {output.narrowPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_FE_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method ppois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method ppois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_PV_bigwig}.bdg.tmp > {output.narrowPeak_PV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_PV_bigwig}.bdg.tmp > {output.narrowPeak_PV_bigwig}.bdg
			bedGraphToBigWig {output.narrowPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_PV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method qpois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method qpois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_qpois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_QV_bigwig}.bdg.tmp > {output.narrowPeak_QV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_QV_bigwig}.bdg.tmp > {output.narrowPeak_QV_bigwig}.bdg
			bedGraphToBigWig {output.narrowPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_QV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			if [ -f {output.narrowPeak_QV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_qpois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_qpois.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.QV.bigwig.bdg
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			fi
		""")


rule Peak_Calling_Broad:
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam",
		processed_index_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{sample, ((?!.*_VS_.*|.*_OVERLAPPED_.*|\\.).)*}.processed.bam.bai",
	output:
		#
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.gz",
		broadPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.bb",
		broadPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.FE.bigwig",
		broadPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.PV.bigwig",
		broadPeak_QV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{sample, ((?!.*_VS_.*|\\.).)*}.broadPeak.QV.bigwig",
		
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {input.processed_bam})
			sample_Name=${{sample_Name%.processed.bam}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling
			#
			if [ ! -f ./Script/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O ./Script/bigBroadPeak.as
			fi
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bam} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_bam} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			printf "%s\\n" "awk '$AWK_COMMAND' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+4 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.broadPeak_bed}.sorted {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted
			awk 'BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+4 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			rm -rf {output.broadPeak_bed}.sorted {output.broadPeak_bed}.tmp
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_FE_bigwig}.bdg.tmp > {output.broadPeak_FE_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_FE_bigwig}.bdg.tmp > {output.broadPeak_FE_bigwig}.bdg
			bedGraphToBigWig {output.broadPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method ppois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method ppois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_PV_bigwig}.bdg.tmp > {output.broadPeak_PV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_PV_bigwig}.bdg.tmp > {output.broadPeak_PV_bigwig}.bdg
			bedGraphToBigWig {output.broadPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method qpois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method qpois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_qpois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_QV_bigwig}.bdg.tmp > {output.broadPeak_QV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_QV_bigwig}.bdg.tmp > {output.broadPeak_QV_bigwig}.bdg
			bedGraphToBigWig {output.broadPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			if [ -f {output.broadPeak_QV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_qpois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.QV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.QV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_qpois.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			fi
		""")


rule Peak_Calling_Broad_Controlled:
	input:
		processed_case_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
	output:
		#
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.gz",
		broadPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.bb",
		broadPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.FE.bigwig",
		broadPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.PV.bigwig",
		broadPeak_QV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.QV.bigwig",
		
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	run:
		shell("""
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {output.broadPeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{GENOME}/{wildcards.design}/peak_calling
			#
			if [ ! -f ./Script/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O ./Script/bigBroadPeak.as
			fi
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_Bam} --control {input.processed_control_Bam} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_case_Bam} --control {input.processed_control_Bam} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			printf "%s\\n" "awk '$AWK_COMMAND' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+4 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "rm -rf {output.broadPeak_bed}.sorted {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted
			awk 'BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+4 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			rm -rf {output.broadPeak_bed}.sorted {output.broadPeak_bed}.tmp
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_FE_bigwig}.bdg.tmp > {output.broadPeak_FE_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_FE_bigwig}.bdg.tmp > {output.broadPeak_FE_bigwig}.bdg
			bedGraphToBigWig {output.broadPeak_FE_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_FE_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method ppois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method ppois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_PV_bigwig}.bdg.tmp > {output.broadPeak_PV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_PV_bigwig}.bdg.tmp > {output.broadPeak_PV_bigwig}.bdg
			bedGraphToBigWig {output.broadPeak_PV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_PV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method qpois" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method qpois
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_qpois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}.bdg.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_QV_bigwig}.bdg.tmp > {output.broadPeak_QV_bigwig}.bdg"  | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}.bdg.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_QV_bigwig}.bdg.tmp > {output.broadPeak_QV_bigwig}.bdg
			bedGraphToBigWig {output.broadPeak_QV_bigwig}.bdg {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_QV_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			if [ -f {output.broadPeak_QV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_qpois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.QV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.QV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_qpois.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.QV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.QV.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak
				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			fi
		""")

