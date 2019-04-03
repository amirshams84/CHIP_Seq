shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Aug-14-2018
# Email: amir.shams84@gmail.com
# Aim: Snakemake recepie for CHIP-Seq Comprehensive analysis
# snakemake --snakefile chip_seq_mouse_new_effort.py --configfile chip_seq_mouse_config.json --debug-dag --cores=50
# snakemake --snakefile chip_seq_mouse.py --configfile chip_seq_mouse_config.json --rulegraph | dot -Tsvg > dag.svg
# ################################### IMPORT ##################################### #


import os
import re
from os.path import join
import sys
import glob
import logging as log
import itertools
import collections
import multiprocessing
# ################################### FUNCTIONS ################################## #


def isFileExist(fname):
	# check the exitence/access/path of a file
	if os.path.isfile(fname) and os.path.exists(fname) and os.access(fname, os.R_OK):
		return True
	else:
		return False


def fix_path(the_Path):
	#
	if the_Path[-1] == "/":
		the_Path = the_Path[:-1]
	return the_Path


def free_memory():
	tot_m, used_m, free_m = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
	return free_m


def get_fq1(wildcards):
	return glob.glob(DATADIR + "/" + wildcards.sample + "*" + config_data_Dict["FORWARD_PATTERN"])


def get_fq2(wildcards):
	return glob.glob(DATADIR + "/" + wildcards.sample + "*" + config_data_Dict["REVERSE_PATTERN"])


def get_processed_tagAlign(wildcards):
	processed_tagAlign_Dict = {}
	processed_tagAlign_Dict["case_processed_tagAlign_List"] = []
	processed_tagAlign_Dict["control_processed_tagAlign_List"] = []
	for Index, study_Dict in config_metadata_Dict.items():
		if study_Dict["Title"] == wildcards.study:
			for case in study_Dict["Case"]:
				processed_tagAlign_Dict["case_processed_tagAlign_List"].append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{case}.processed.tagAlign.gz".format(study=wildcards.study, case=case))
			for control in study_Dict["Control"]:
				processed_tagAlign_Dict["control_processed_tagAlign_List"].append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{control}.processed.tagAlign.gz".format(study=wildcards.study, control=control))
	return processed_tagAlign_Dict


def get_case_processed_tagAlign(wildcards):
	processed_tagAlign_List = []
	for Index, study_Dict in config_metadata_Dict.items():
		if study_Dict["Title"] == wildcards.study:
			for case in study_Dict["Case"]:
				processed_tagAlign_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample}.processed.tagAlign.gz".format(study=wildcards.study, sample=case))
	return processed_tagAlign_List


def get_control_processed_tagAlign(wildcards):
	processed_tagAlign_List = []
	for Index, study_Dict in config_metadata_Dict.items():
		if study_Dict["Title"] == wildcards.study:
			for control in study_Dict["Control"]:
				processed_tagAlign_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample}.processed.tagAlign.gz".format(study=wildcards.study, sample=control))
	return processed_tagAlign_List

#def get_tagAlign(wildcards):
	#processed_tagAlign_List = []
	#processed_tagAlign_List.extend()

# ################################### CONFIGURATION ############################## #


localrules: all
configfile: "chip_seq_mouse_config.json"
config_general_Dict = config["GENERAL"]
config_data_Dict = config["DATA"]
config_metadata_Dict = config["METADATA"]
config_pre_processing_Dict = config["PRE_PROCESSING"]
config_alignment_Dict = config["ALIGNMENT"]
config_post_alignment_Dict = config["POST_ALIGNMENT"]
config_peak_calling_Dict = config["PEAK_CALLING"]
config_reference_Dict = config["REFERENCE"]
config_peak_analysis_Dict = config["PEAK_ANALYSIS"]
#
WORKDIR = fix_path(config_general_Dict["WORK_DIR"])
DATADIR = fix_path(config_general_Dict["DATA_DIR"])
EXECDIR = fix_path(config_general_Dict["EXEC_DIR"])
REFDIR = fix_path(config_general_Dict["REF_DIR"])
#
TITLE = config_general_Dict["TITLE"]
PROCESSORS = multiprocessing.cpu_count()
#MEMORY = free_memory()
MEMORY = 4000
#
CONDA_INIT = "source /data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/.conda_env"
CONDA_ENV_PY2 = "/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2"
CONDA_ENV_PY3 = "/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py3"
ACTIVATE_ENV_PY2 = "source activate /data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2"
ACTIVATE_ENV_PY3 = "source activate /data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py3"
# ################################### WILDCARDS ################################## #

case_List = []
case_PRE_PROCESSING_List = []
case_ALIGNMENT_List = []
case_POST_ALIGNMENT_List = []
case_NARROW_PEAK_CALLING_List = []
control_List = []
control_PRE_PROCESSING_List = []
control_ALIGNMENT_List = []
control_POST_ALIGNMENT_List = []
control_NARROW_PEAK_CALLING_List = []
POOLING_REPLICATE_List = []
NARROW_PEAK_CALLING_POOLED_List = []
NARROW_PEAK_CALLING_POOLED_VERSUS_List = []
for Index, study_Dict in config_metadata_Dict.items():
	#
	POOLING_REPLICATE_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{pooled}.processed.tagAlign.gz".format(study=study_Dict["Title"], pooled="-POOLED_CASE-".join(study_Dict["Case"])))
	POOLING_REPLICATE_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{pooled}.processed.tagAlign.gz".format(study=study_Dict["Title"], pooled="-POOLED_CONTROL-".join(study_Dict["Control"])))
	NARROW_PEAK_CALLING_POOLED_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{pooled}.narrowPeak.gz".format(study=study_Dict["Title"], pooled="-POOLED_CASE-".join(study_Dict["Case"])))
	NARROW_PEAK_CALLING_POOLED_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{pooled}.narrowPeak.gz".format(study=study_Dict["Title"], pooled="-POOLED_CONTROL-".join(study_Dict["Control"])))
	NARROW_PEAK_CALLING_POOLED_VERSUS_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{pooled_case}_VS_{pooled_control}.narrowPeak.gz".format(study=study_Dict["Title"], pooled_case="-POOLED_CASE-".join(study_Dict["Case"]), pooled_control="-POOLED_CONTROL-".join(study_Dict["Control"])))
	for case in study_Dict["Case"]:
		#
		if case not in case_List:
			case_List.append("{case}".format(case=case))
			case_PRE_PROCESSING_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{case}_R1.processed.fastq.gz".format(study=study_Dict["Title"], case=case))
			case_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{case}.bam".format(study=study_Dict["Title"], case=case))
			case_POST_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{case}.processed.tagAlign.gz".format(study=study_Dict["Title"], case=case))
			case_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}.narrowPeak.gz".format(study=study_Dict["Title"], case=case))
	for control in study_Dict["Control"]:
		#
		if control not in control_List:
			control_List.append("{control}".format(control=control))
			control_PRE_PROCESSING_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{control}_R1.processed.fastq.gz".format(study=study_Dict["Title"], control=control))
			control_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{control}.bam".format(study=study_Dict["Title"], control=control))
			control_POST_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{control}.processed.tagAlign.gz".format(study=study_Dict["Title"], control=control))
			control_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{control}.narrowPeak.gz".format(study=study_Dict["Title"], control=control))


sample_List = case_List + control_List
PRE_PROCESSING_List = case_PRE_PROCESSING_List + control_PRE_PROCESSING_List
ALIGNMENT_List = case_ALIGNMENT_List + control_ALIGNMENT_List
POST_ALIGNMENT_List = case_POST_ALIGNMENT_List + control_POST_ALIGNMENT_List
NARROW_PEAK_CALLING_List = case_NARROW_PEAK_CALLING_List + control_NARROW_PEAK_CALLING_List


case_vs_control_NARROW_PEAK_CALLING_List = []
for Index, study_Dict in config_metadata_Dict.items():
	#
	for case in study_Dict["Case"]:
		#
		for control in study_Dict["Control"]:
			#
			case_vs_control_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.narrowPeak.gz".format(study=study_Dict["Title"], case=case, control=control))

#wildcard_constraints:
#	sample = "[^/.(?!_VS_)]+",
#	case = "[^/.]+",
#	control = "[^/.]+",

# ################################### RULES ###################################### #
rule all:
	input: PRE_PROCESSING_List + ALIGNMENT_List + POST_ALIGNMENT_List + NARROW_PEAK_CALLING_List + case_vs_control_NARROW_PEAK_CALLING_List + POOLING_REPLICATE_List + NARROW_PEAK_CALLING_POOLED_List + NARROW_PEAK_CALLING_POOLED_VERSUS_List


rule PRE_PROCESSING:
	input:
		fq1 = get_fq1,
		fq2 = get_fq2,
	output:
		processed_R1 = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}_R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}_R2.processed.fastq.gz",
		discarded_R1 = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}_R1.discarded.fastq.gz",
		discarded_R2 = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}_R2.discarded.fastq.gz",
		cutadapt_report = WORKDIR + "/" + TITLE + "/{study}/{replicate}/QUALITY_CONTROL/{sample, [^/.(?!_VS_)]+}.cutadapt_report.txt",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}.PRE_PROCESSING.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}.PRE_PROCESSING.stderr",
	priority: 99
	threads: PROCESSORS
	message: "PRE_PROCESSING"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			printf "#############################################################################\\n"
			printf "PRE_PROCESSING\\n"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "Step1: REMOVING SHORT,LOW_COMPLEXITY READS and TRIMMING THE ADAPTER\\n"
			printf "\\tAPPLICATION1: cutadapt version:%s\\n" "$(cutadapt --version)"
			printf "\\tPARAMETERS: %s\\n" "{config_pre_processing_Dict[CUTADAPT]}"
			printf "\\tINPUT1: %s\\n" "{input.fq1}"
			printf "\\tINPUT2: %s\\n" "{input.fq2}"
			printf "\\tOUTPUT1: %s\\n" "{output.processed_R1}"
			printf "\\tOUTPUT2: %s\\n" "{output.processed_R2}"
			printf "\\tOUTPUT3: %s\\n" "{output.discarded_R1}"
			printf "\\tOUTPUT4: %s\\n" "{output.discarded_R2}"
			printf "\\tOUTPUT5: %s\\n" "{output.cutadapt_report}"
			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "cutadapt {config_pre_processing_Dict[CUTADAPT]} --too-short-output={output.discarded_R1} --too-short-paired-output={output.discarded_R2} --output={output.processed_R1} --paired-output={output.processed_R2} {input.fq1} {input.fq2}"
			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			cutadapt {config_pre_processing_Dict[CUTADAPT]} --too-short-output={output.discarded_R1} --too-short-paired-output={output.discarded_R2} --output={output.processed_R1} --paired-output={output.processed_R2} \
			{input.fq1} {input.fq2} > {output.cutadapt_report} 2>> {log.stderr}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "#############################################################################\\n"
			""")


rule ALIGNMENT:
	input:
		processed_R1 = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}_R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PRE_PROCESSING/{sample, [^/.(?!_VS_)]+}_R2.processed.fastq.gz"
	output:
		alignment_Bam = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.bam",
		alignment_Report = WORKDIR + "/" + TITLE + "/{study}/{replicate}/QUALITY_CONTROL/{sample, [^/.(?!_VS_)]+}.alignment_report.txt",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.ALIGNMENT.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.ALIGNMENT.stderr",

	priority: 98
	threads: PROCESSORS
	message: "ALIGNMENT"
	resources:
		mem_mb = MEMORY

	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			R1=$(basename {input.processed_R1})
			printf "#############################################################################\\n"
			printf "ALIGNMENT\\n"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "Step1: ALIGNMENT AGAINST {config_alignment_Dict[TYPE]} GENOME\\n"
			printf "\\tAPPLICATION1: bowtie2 version:%s\\n" "$(bowtie2 --version | head -1 | cut -f3 -d' ')"
			printf "\\tPARAMETERS: %s\\n" "{config_alignment_Dict[BOWTIE2]}"
			printf "\\tREFERENCE: %s\\n" "{config_reference_Dict[MM10][Bowtie2Index]}"
			printf "\\tINPUT1: %s\\n" "{input.processed_R1}"
			printf "\\tINPUT2: %s\\n" "{input.processed_R2}"
			printf "\\tOUTPUT1: %s\\n" "{output.alignment_Bam}"
			printf "\\tOUTPUT2: %s\\n" "{output.alignment_Report}"
			printf "\\tOUTPUT3: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R1.processed.unmapped.fastq.gz"
			printf "\\tOUTPUT4: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R2.processed.unmapped.fastq.gz"
			printf "\\tOUTPUT5: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R1.processed.mapped.fastq.gz"
			printf "\\tOUTPUT6: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R2.processed.mapped.fastq.gz"
			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "bowtie2 {config_alignment_Dict[BOWTIE2]} --threads {threads} -x {config_reference_Dict[MM10][Bowtie2Index]} -1 {input.processed_R1} -2 {input.processed_R2} --un-conc-gz {WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R%.processed.unmapped.fastq.gz --al-conc-gz {WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R%.processed.mapped.fastq.gz 2> {output.alignment_Report} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{TITLE}/{wildcards.study}/ALIGNMENT/{wildcards.sample} -o {output.alignment_Bam} -"
			printf "\\t%s\\n" "samtools index -@ {threads} -b {output.alignment_Bam}"
			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			bowtie2 {config_alignment_Dict[BOWTIE2]} --threads {threads} -x {config_reference_Dict[MM10][Bowtie2Index]} -1 {input.processed_R1} -2 {input.processed_R2} \
			--un-conc-gz {WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R%.processed.unmapped.fastq.gz \
			--al-conc-gz {WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/{wildcards.sample}_R%.processed.mapped.fastq.gz \
			2> {output.alignment_Report} | samtools sort --threads {threads} -O bam -T {WORKDIR}/{TITLE}/{wildcards.study}/ALIGNMENT/{wildcards.sample} -o {output.alignment_Bam} - >> {log.stdout}
			samtools index -@ 1 -b {output.alignment_Bam} >> {log.stdout} 2>> {log.stderr}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "#############################################################################\\n"
			""")


rule POST_ALIGNMENT:
	input:
		alignment_Bam = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.bam",
	output:
		#PICARD
		marked_duplicate_Bam = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.marked_duplicate.bam",
		picard_report = WORKDIR + "/" + TITLE + "/{study}/{replicate}/QUALITY_CONTROL/{sample, [^/.(?!_VS_)]+}.picard_report.txt",
		flagstat = WORKDIR + "/" + TITLE + "/{study}/{replicate}/QUALITY_CONTROL/{sample, [^/.(?!_VS_)]+}.flagstat.txt",
		#SAMTOOLS
		processed_Bam = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.processed.bam",
		processed_Flagstat = WORKDIR + "/" + TITLE + "/{study}/{replicate}/QUALITY_CONTROL/{sample, [^/.(?!_VS_)]+}.processed.flagstat.txt",
		nmsort_Bam = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.processed.nmsort.bam",
		Bed = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.bed.gz",
		processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.processed.tagAlign.gz",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.POST_ALIGNMENT.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.POST_ALIGNMENT.stderr",
	priority: 97
	threads: PROCESSORS
	message: "POST_ALIGNMENT"
	resources:
		mem_mb = MEMORY

	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			Sample_Name=$(basename {input.alignment_Bam})
			Sample_Name=${{Sample_Name%.bam}}
			printf "#############################################################################\\n"
			printf "POST_ALIGNMENT\\n"
			printf "++++++++++++++++++++++++++++++++++++++++\\n"
			printf "Step1: MARKING DUPLICATE\\n"
			printf "\\tAPPLICATION1: picard version:%s\\n" "$(picard MarkDuplicates --version 2>&1 | head -n 1)"
			printf "\\tPARAMETERS: %s\\n" "{config_post_alignment_Dict[PICARD]}"
			printf "\\tINPUT1: %s\\n" "{input.alignment_Bam}"
			printf "\\tOUTPUT1: %s\\n" "{output.marked_duplicate_Bam}"
			printf "\\tOUTPUT2: %s\\n" "{output.picard_report}"
			export _JAVA_OPTIONS="-Xms{MEMORY}M -Xmx{MEMORY}M -XX:ParallelGCThreads={threads}"
			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "picard MarkDuplicates INPUT={input.alignment_Bam} OUTPUT={output.marked_duplicate_Bam} METRICS_FILE={output.picard_report} {config_post_alignment_Dict[PICARD]}"
			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			picard MarkDuplicates INPUT={input.alignment_Bam} OUTPUT={output.marked_duplicate_Bam} METRICS_FILE={output.picard_report} {config_post_alignment_Dict[PICARD]} >> {log.stdout} 2>> {log.stderr}
			samtools flagstat --threads {threads} {output.marked_duplicate_Bam} > {output.flagstat}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "Step2: FILTERING UNWANTED MAPPING\\n"
			printf "\\tAPPLICATION1: samtools version:%s\\n" "$(samtools --version | head -1 | cut -f2 -d' ')"
			printf "\\tPARAMETERS: %s\\n" "{config_post_alignment_Dict[SAMTOOLS]}"
			printf "\\tINPUT1: %s\\n" "{output.marked_duplicate_Bam}"
			printf "\\tOUTPUT1: %s\\n" "{output.processed_Bam}"
			printf "\\tCOMMAND_LINE:\\n"
			AWK_COMMAND='{{if($3 !~ /chrUn/ && $3 != "ChrEBV" && $3 !~ /random/ ){{print $0}}}}'
			printf "\\t%s | %s | %s\\n" "samtools view --threads {threads} -h {output.marked_duplicate_Bam}" "awk '$AWK_COMMAND'" "samtools view --threads {threads} -F 1804 -f 2 -q 30 -b - -o {output.processed_Bam}"
			printf "\\t%s\\n" "samtools index -@ {threads} -b {output.processed_Bam}"
			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			samtools view --threads {threads} -h {output.marked_duplicate_Bam} | awk '{{if($3 !~ /chrUn/ && $3 != "ChrEBV" && $3 !~ /random/ ){{print $0}}}}' | \
			samtools view --threads {threads} -F 1804 -f 2 -q 30 -b - -o {output.processed_Bam} >> {log.stdout} 2>> {log.stderr}
			samtools index -@ {threads} -b {output.processed_Bam} >> {log.stdout} 2>> {log.stderr}
			samtools flagstat --threads {threads} {output.processed_Bam} > {output.processed_Flagstat}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step3: TAGALIGN FORMATTING and ADJUST ALIGNMENT"
			printf "\\tAPPLICATION1: samtools version:%s\\n" "$(samtools --version | head -1 | cut -f2 -d' ')"
			printf "\\tAPPLICATION2: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tINPUT1: %s\\n" "{output.processed_Bam}"
			printf "\\tOUTPUT1: %s\\n" "{output.nmsort_Bam}"
			printf "\\tOUTPUT2: %s\\n" "{output.Bed}"
			printf "\\tOUTPUT3: %s\\n" "{output.processed_tagAlign}"
			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "samtools sort --threads {threads} -n {output.processed_Bam} -o {output.nmsort_Bam}"
			printf "\\t%s\\n" "bedtools bamtobed -bedpe -mate1 -i {output.nmsort_Bam} | gzip -nc > {output.Bed}"
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",$1,$2,$3,$9,$4,$5,$6,$10}}'
			printf "\\t%s | %s | %s\\n" "zcat {output.Bed}" "awk '$AWK_COMMAND'" "gzip -nc > {output.processed_tagAlign}"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			samtools sort --threads {threads} -n {output.processed_Bam} -o {output.nmsort_Bam} >> {log.stdout} 2>> {log.stderr}
			bedtools bamtobed -bedpe -mate1 -i {output.nmsort_Bam} | gzip -nc > {output.Bed}
			zcat {output.Bed} | awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",$1,$2,$3,$9,$4,$5,$6,$10}}' | gzip -nc > {output.processed_tagAlign}

			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			echo "#############################################################################"
			""")

rule POOLING_REPLICATE:
	input:
		case_processed_tagAlign_List = expand(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample}.processed.tagAlign.gz", study=study_Dict["Title"], sample=case_List),
		control_processed_tagAlign_List = expand(WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample}.processed.tagAlign.gz", study=study_Dict["Title"], sample=control_List)
	output:
		pooled_case_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{pooled, .*-POOLED_CASE-.*}.processed.tagAlign.gz",
		pooled_control_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{pooled, .*-POOLED_CONTROL-.*}.processed.tagAlign.gz",
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			echo "pooling cases"
			printf "\\tINPUT1: [%s]\\n" "{input.case_processed_tagAlign_List}"
			zcat {input.case_processed_tagAlign_List} | gzip -nc > {output.pooled_case_processed_tagAlign}
			zcat {input.control_processed_tagAlign_List} | gzip -nc > {output.pooled_control_processed_tagAlign}
			""")


rule NARROW_PEAK_CALLING:
	input:
		processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, [^/.(?!_VS_)]+}.processed.tagAlign.gz",
	output:
		narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{sample, [^/.(?!_VS_)]+}.narrowPeak.gz",
	priority: 95
	message: "NARROW_PEAK_CALLING"
	threads: PROCESSORS

	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			sample_Name=$(basename {input.processed_tagAlign})
			sample_Name=${{sample_Name%.processed.tagAlign.gz}}
			echo "macs2 callpeak --treatment {input.processed_tagAlign} -f BED --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			macs2 callpeak --treatment {input.processed_tagAlign} -f BED --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak | awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz | sed '/^\\(chr\\)/!d' | sort -k1,1V -k2,2n > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.hammock
			/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2/script/narrowpeak2hammock.py {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.hammock

			bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )  | grep -P 'chr[\\dXY]+[\\t]' | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz | awk 'BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}' | sort -grk8 | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz | sed '/^\\(chr\\)/!d' | sort -k1,1V -k2,2n > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.hammock
			/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2/script/narrowpeak2hammock.py {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.hammock

			sval=$(wc -l <(zcat -f {input.processed_tagAlign}) | awk '{{printf "%f", $1/1000000}}')
			echo "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois -S ${{sval}} "
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name \
			--outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois -S ${{sval}}

			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.sorted.bedgraph
			mv {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.sorted.bedgraph {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph.gz

			bedGraphToBigWig {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bigwig
			""")

rule NARROW_PEAK_CALLING_POOLED:
	input:
		processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, (?!.*_VS_).*-POOLED_.*}.processed.tagAlign.gz",
	output:
		narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{sample, (?!.*_VS_).*-POOLED_.*}.narrowPeak.gz",
	priority: 95
	message: "NARROW_PEAK_CALLING_POOLED"
	threads: PROCESSORS

	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			sample_Name=$(basename {input.processed_tagAlign})
			sample_Name=${{sample_Name%.processed.tagAlign.gz}}
			echo "macs2 callpeak --treatment {input.processed_tagAlign} -f BED --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			macs2 callpeak --treatment {input.processed_tagAlign} -f BED --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak | awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz | sed '/^\\(chr\\)/!d' | sort -k1,1V -k2,2n > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.hammock
			/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2/script/narrowpeak2hammock.py {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.hammock

			bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )  | grep -P 'chr[\\dXY]+[\\t]' | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz | awk 'BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}' | sort -grk8 | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz | sed '/^\\(chr\\)/!d' | sort -k1,1V -k2,2n > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.hammock
			/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2/script/narrowpeak2hammock.py {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.hammock

			sval=$(wc -l <(zcat -f {input.processed_tagAlign}) | awk '{{printf "%f", $1/1000000}}')
			echo "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois -S ${{sval}} "
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name \
			--outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois -S ${{sval}}

			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.sorted.bedgraph
			mv {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.sorted.bedgraph {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph.gz

			bedGraphToBigWig {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bigwig
			""")


rule NARROW_PEAK_CALLING_VERSUS:
	input:
		case_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{case, [^/.(?!_VS_)]+(.*-POOLED_.*)*}.processed.tagAlign.gz",
		control_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{control, [^/.(?!_VS_)]+(.*-POOLED_.*)*}.processed.tagAlign.gz",
	output:
		versus_narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case, [^/.(?!_VS_)]+(.*-POOLED_.*)*}_VS_{control, [^/.(?!_VS_)]+(.*-POOLED_.*)*}.narrowPeak.gz",
	priority: 95
	message: "NARROW_PEAK_CALLING_VERSUS"
	threads: PROCESSORS

	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			sample_Name=$(basename {output.versus_narrowPeak})
			sample_Name=${{sample_Name%.narrowPeak.gz}}
			echo "macs2 callpeak --treatment {input.case_processed_tagAlign} --control {input.control_processed_tagAlign}  -f BED --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			macs2 callpeak --treatment {input.case_processed_tagAlign} --control {input.control_processed_tagAlign}  -f BED --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_peaks.narrowPeak | awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz | sed '/^\\(chr\\)/!d' | sort -k1,1V -k2,2n > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.hammock
			/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2/script/narrowpeak2hammock.py {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.hammock

			bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.narrowPeak.gz ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )  | grep -P 'chr[\\dXY]+[\\t]' | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.blk_filt.narrowPeak.gz | awk 'BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}' | sort -grk8 | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.gz | sed '/^\\(chr\\)/!d' | sort -k1,1V -k2,2n > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.hammock
			/data/shamsaddinisha/ATAC_Seq/WORKDIR/EXECDIR/Miniconda3/envs/atac_seq_py2/script/narrowpeak2hammock.py {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}.processed.narrowPeak.hammock

			sval=$(wc -l <(zcat -f {input.case_processed_tagAlign}) | awk '{{printf "%f", $1/1000000}}')
			echo "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois -S ${{sval}} "
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{sample_Name}}_control_lambda.bdg --o-prefix $sample_Name \
			--outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois -S ${{sval}}

			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph
			LC_COLLATE=C sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.sorted.bedgraph
			mv {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.sorted.bedgraph {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph
			bgzip -c {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph.gz

			bedGraphToBigWig {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bedgraph {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{sample_Name}}.bigwig
			""")

'''

rule NARROW_PEAK_CALLING_POOLED:
	input:
		processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{sample, .*-POOLED_.*}.processed.tagAlign.gz",
	output:
		narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{sample, .*-POOLED_.*}.narrowPeak.gz",
	priority: 95
	message: "NARROW_PEAK_CALLING"
	threads: PROCESSORS

	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			echo {input.processed_tagAlign}
			touch {output.narrowPeak}
			""")

rule NARROW_PEAK_CALLING_POOLED_VERSUS:
	input:
		pooled_case_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{pooled_case, .*-POOLED_CASE-.*}.processed.tagAlign.gz",
		pooled_control_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{pooled_control, .*-POOLED_CONTROL-.*}.processed.tagAlign.gz",
	output:
		narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{pooled_case, .*-POOLED_CONTROL-.*}_VS_{pooled_control, .*-POOLED_CONTROL-.*}.narrowPeak.gz",
	priority: 95
	message: "NARROW_PEAK_CALLING"
	threads: PROCESSORS

	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			echo {input.pooled_case_processed_tagAlign}
			echo {input.pooled_control_processed_tagAlign}
			touch {output.narrowPeak}
			""")


rule NARROW_PEAK_CALLING:
	input:
		case_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{case}.processed.tagAlign.gz",
		control_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/{replicate}/ALIGNMENT/{control}.processed.tagAlign.gz",
	output:
		case_vs_control_narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.narrowPeak.gz",

	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.NARROW_PEAK_CALLING.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.NARROW_PEAK_CALLING.stderr",
	priority: 95
	message: "NARROW_PEAK_CALLING"
	threads: PROCESSORS

	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			case_Name=$(basename {input.case_processed_tagAlign})
			case_Name=${{case_Name%.processed.tagAlign.gz}}
			printf "%s\\n" "#############################################################################"
			printf "%s\\n" "PEAK CALLING"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step1: NARROW PEAK CALLING OF $case_Name"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tPARAMETERS: %s\\n" "{config_peak_calling_Dict[MACS2_NARROW]}"
			printf "\\tINPUT1: %s\\n" "{input.case_processed_tagAlign}"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/PEAK_CALLING/MACS2/${{case_Name}}.narrowPeak.gz"
			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "macs2 callpeak --treatment {input.case_processed_tagAlign} -f BED --name $case_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}'
			printf "\\t%s | %s | %s\\n" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_peaks.narrowPeak" "awk '$AWK_COMMAND'" "gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.narrowPeak.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.case_processed_tagAlign} -f BED --name $case_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/ >> {log.stdout} 2>> {log.stderr}
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_peaks.narrowPeak | awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.narrowPeak.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			control_Name=$(basename {input.control_processed_tagAlign})
			control_Name=${{control_Name%.processed.tagAlign.gz}}
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step2: NARROW PEAK CALLING OF $control_Name"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tPARAMETERS: %s\\n" "{config_peak_calling_Dict[MACS2_NARROW]}"
			printf "\\tINPUT1: %s\\n" "{input.control_processed_tagAlign}"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/PEAK_CALLING/MACS2/${{control_Name}}.narrowPeak.gz"
			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "macs2 callpeak --treatment {input.control_processed_tagAlign} -f BED --name $control_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}'
			printf "\\t%s | %s | %s\\n" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_peaks.narrowPeak" "awk '$AWK_COMMAND'" "gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.narrowPeak.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.control_processed_tagAlign} -f BED --name $control_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/ >> {log.stdout} 2>> {log.stderr}
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_peaks.narrowPeak | awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.narrowPeak.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			case_vs_control_Name=$(basename {output.case_vs_control_narrowPeak})
			case_vs_control_Name=${{case_vs_control_Name%.narrowPeak.gz}}
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step3: NARROW PEAK CALLING OF ${{case_Name}}_VS_${{control_Name}}"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tPARAMETERS: %s\\n" "{config_peak_calling_Dict[MACS2_NARROW]}"
			printf "\\tINPUT1: %s\\n" "{input.case_processed_tagAlign}"
			printf "\\tINPUT2: %s\\n" "{input.control_processed_tagAlign}"
			printf "\\tOUTPUT1: %s\\n" "{output.case_vs_control_narrowPeak}"
			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "macs2 callpeak --treatment {input.case_processed_tagAlign} --control {input.control_processed_tagAlign} -f BED --name $case_vs_control_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}'
			printf "\\t%s | %s | %s\\n" "sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_peaks.narrowPeak" "awk '$AWK_COMMAND'" "gzip -nc > {output.case_vs_control_narrowPeak}"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 callpeak --treatment {input.case_processed_tagAlign} --control {input.control_processed_tagAlign} -f BED --name $case_vs_control_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/ >> {log.stdout} 2>> {log.stderr}
			sort -k 8gr,8gr {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_peaks.narrowPeak | awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}' | gzip -nc > {output.case_vs_control_narrowPeak}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "#############################################################################"
			""")


rule POST_PEAK_CALLING:
	input:
		narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.narrowPeak.gz",
	output:
		blk_filt_narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.blk_filt.narrowPeak.gz",
		processed_narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.processed.narrowPeak.gz",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.POST_PEAK_CALLING.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.POST_PEAK_CALLING.stderr",
	priority: 95
	threads: PROCESSORS
	message: "POST_PEAK_CALLING"
	resources:
		mem_mb = MEMORY

	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			mkdir -p {REFDIR}
			if [ ! -f {REFDIR}/mm10.chrom.sizes ]; then
				wget {config_reference_Dict[MM10][CHROMOSOME_SIZE]} -O {REFDIR}/mm10.chrom.sizes >> {log.stdout} 2>> {log.stderr}
			fi

			if [ ! -f {REFDIR}/mm10.blacklist.bed.gz ]; then
				wget {config_reference_Dict[MM10][BLACKLIST]} -O {REFDIR}/mm10.blacklist.bed.gz >> {log.stdout} 2>> {log.stderr}
			fi
			case_vs_control_Name=$(basename {input.narrowPeak})
			case_Name=${{case_vs_control_Name/%_VS*/}}
			printf "%s\\n" "#############################################################################"
			printf "%s\\n" "POST PEAK CALLING: $case_Name"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step1: FILTERING BLACKLIST REGIONS"
			printf "\\tAPPLICATION1: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.narrowPeak.gz"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.blk_filt.narrowPeak.gz"
			printf "\\tCOMMAND_LINE:\\n"
			AWK_COMMAND='BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			printf "\\t%s | %s | %s | %s\\n" "bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.narrowPeak.gz ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )"  "grep -P 'chr[\\dXY]+[\\t]'" "awk '$AWK_COMMAND'" "gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.blk_filt.narrowPeak.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.narrowPeak.gz ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )  | grep -P 'chr[\\dXY]+[\\t]' | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.blk_filt.narrowPeak.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step2: FILTERING SIGNIFICANT PEAKS"
			printf "\\tAPPLICATION1: zcat version:%s\\n" "$(zcat --version | head -1 | cut -f3 -d' ')"
			printf "\\tPARAMETERS: P-VALUE_THRESHOLD=%s\\n" "{config_peak_calling_Dict[PVALUE]}"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.blk_filt.narrowPeak.gz"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.processed.narrowPeak.gz"
			printf "\\tCOMMAND_LINE:\\n"
			AWK_COMMAND='BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}'
			printf "\\t%s | %s | %s | %s\\n" "zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.blk_filt.narrowPeak.gz" "awk '$AWK_COMMAND'" "sort -grk8" "gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.processed.narrowPeak.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.blk_filt.narrowPeak.gz | awk 'BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}' | sort -grk8 | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}.processed.narrowPeak.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "#############################################################################"


			case_vs_control_Name=$(basename {input.narrowPeak})
			control_Name=$(echo "${{case_vs_control_Name}}" | sed -e "s/^${{case_Name}}_VS_//" -e "s/.narrowPeak.gz$//")
			printf "%s\\n" "#############################################################################"
			printf "%s\\n" "POST PEAK CALLING: $control_Name"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step1: FILTERING BLACKLIST REGIONS"
			printf "\\tAPPLICATION1: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.narrowPeak.gz"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.blk_filt.narrowPeak.gz"
			printf "\\tCOMMAND_LINE:\\n"
			AWK_COMMAND='BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			printf "\\t%s | %s | %s | %s\\n" "bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.narrowPeak.gz ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )"  "grep -P 'chr[\\dXY]+[\\t]'" "awk '$AWK_COMMAND'" "gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.blk_filt.narrowPeak.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			bedtools intersect -v -a <(zcat -f {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.narrowPeak.gz ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )  | grep -P 'chr[\\dXY]+[\\t]' | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.blk_filt.narrowPeak.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step2: FILTERING SIGNIFICANT PEAKS"
			printf "\\tAPPLICATION1: zcat version:%s\\n" "$(zcat --version | head -1 | cut -f3 -d' ')"
			printf "\\tPARAMETERS: P-VALUE_THRESHOLD=%s\\n" "{config_peak_calling_Dict[PVALUE]}"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.blk_filt.narrowPeak.gz"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.processed.narrowPeak.gz"
			printf "\\tCOMMAND_LINE:\\n"
			AWK_COMMAND='BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}'
			printf "\\t%s | %s | %s | %s\\n" "zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.blk_filt.narrowPeak.gz" "awk '$AWK_COMMAND'" "sort -grk8" "gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.processed.narrowPeak.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			zcat {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.blk_filt.narrowPeak.gz | awk 'BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}' | sort -grk8 | gzip -nc > {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}.processed.narrowPeak.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "#############################################################################"
			printf "%s\\n" "POST PEAK CALLING: ${{case_Name}}_VS_${{control_Name}}"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step1: FILTERING BLACKLIST REGIONS"
			printf "\\tAPPLICATION1: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tINPUT1: %s\\n" "{input.narrowPeak}"
			printf "\\tOUTPUT1: %s\\n" "{output.blk_filt_narrowPeak}"
			printf "\\tCOMMAND_LINE:\\n"
			AWK_COMMAND='BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			printf "\\t%s | %s | %s | %s\\n" "bedtools intersect -v -a <(zcat -f {input.narrowPeak} ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )"  "grep -P 'chr[\\dXY]+[\\t]'" "awk '$AWK_COMMAND'" "gzip -nc > {output.blk_filt_narrowPeak}"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			bedtools intersect -v -a <(zcat -f {input.narrowPeak} ) -b <(zcat -f {REFDIR}/mm10.blacklist.bed.gz )  | grep -P 'chr[\\dXY]+[\\t]' | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' | gzip -nc > {output.blk_filt_narrowPeak}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step2: FILTERING SIGNIFICANT PEAKS"
			printf "\\tAPPLICATION1: zcat version:%s\\n" "$(zcat --version | head -1 | cut -f3 -d' ')"
			printf "\\tPARAMETERS: P-VALUE_THRESHOLD=%s\\n" "{config_peak_calling_Dict[PVALUE]}"
			printf "\\tINPUT1: %s\\n" "{output.blk_filt_narrowPeak}"
			printf "\\tOUTPUT1: %s\\n" "{output.processed_narrowPeak}"
			printf "\\tCOMMAND_LINE:\\n"
			AWK_COMMAND='BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}'
			printf "\\t%s | %s | %s | %s\\n" "zcat {output.blk_filt_narrowPeak}" "awk '$AWK_COMMAND'" "sort -grk8" "gzip -nc > {output.processed_narrowPeak}"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			zcat {output.blk_filt_narrowPeak} | awk 'BEGIN{{FS="\\t"}} {{if ($8>-log({config_peak_calling_Dict[PVALUE]})/log(10)) print}}' | sort -grk8 | gzip -nc > {output.processed_narrowPeak}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "#############################################################################"
			""")


rule SIGNAL_PROCESSING:
	input:
		#NARROWPEAK
		narrowPeak = WORKDIR + "/" + TITLE + "/{study}/{replicate}/PEAK_CALLING/MACS2/{case}_VS_{control}.narrowPeak.gz",
	output:
		#fold_enrichment SIGNAL
		fold_enrichment_signal_Bedgraph = WORKDIR + "/" + TITLE + "/{study}/{replicate}/SIGNAL_PROCESSING/MACS2/{case}_VS_{control}.fold_enrichment_signal.bedgraph",
		fold_enrichment_sorted_signal_Bedgraph_gz = WORKDIR + "/" + TITLE + "/{study}/{replicate}/SIGNAL_PROCESSING/MACS2/{case}_VS_{control}.fold_enrichment_signal.bedgraph.gz",
		#Poisson Pvalue SIGNAL
		poisson_pvalue_signal_Bedgraph = WORKDIR + "/" + TITLE + "/{study}/{replicate}/SIGNAL_PROCESSING/MACS2/{case}_VS_{control}.poisson_pvalue_signal.bedgraph",
		poisson_pvalue_sorted_signal_Bedgraph_gz = WORKDIR + "/" + TITLE + "/{study}/{replicate}/SIGNAL_PROCESSING/MACS2/{case}_VS_{control}.poisson_pvalue_signal.bedgraph.gz",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/{replicate}/SIGNAL_PROCESSING/MACS2/{case}_VS_{control}.SIGNAL_PROCESSING.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/{replicate}/SIGNAL_PROCESSING/MACS2/{case}_VS_{control}.SIGNAL_PROCESSING.stderr",
	priority: 95
	threads: PROCESSORS
	message: "SIGNAL_PROCESSING"
	resources:
		mem_mb = MEMORY

	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			mkdir -p {REFDIR}
			if [ ! -f {REFDIR}/mm10.chrom.sizes ]; then
				wget {config_reference_Dict[MM10][CHROMOSOME_SIZE]} -O {REFDIR}/mm10.chrom.sizes >> {log.stdout} 2>> {log.stderr}
			fi

			if [ ! -f {REFDIR}/mm10.blacklist.bed.gz ]; then
				wget {config_reference_Dict[MM10][BLACKLIST]} -O {REFDIR}/mm10.blacklist.bed.gz >> {log.stdout} 2>> {log.stderr}
			fi
			case_vs_control_Name=$(basename {input.narrowPeak})
			case_Name=${{case_vs_control_Name/%_VS*/}}
			printf "%s\\n" "#############################################################################"
			printf "%s\\n" "SIGNAL PROCESSING: ${{case_Name}}"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step1: BUILDING FOLD ENRICHMENT SIGNAL TRACKS"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tAPPLICATION2: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tAPPLICATION3: ucsc-bedclip version:%s\\n" "366"
			printf "\\tAPPLICATION4: ucsc-bedgraphtobigwig version:%s\\n" "366"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_treat_pileup.bdg"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph"
			printf "\\tOUTPUT2: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph.gz"

			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_control_lambda.bdg --o-prefix $case_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m FE "
			printf "\\t%s | %s\\n" "slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}_FE.bdg -g {REFDIR}/mm10.chrom.sizes -b 0" "bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph"
			printf "\\t%s\\n" "sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph | bgzip -c >  {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph.gz"
			printf "\\t%s\\n" "tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_control_lambda.bdg \
			--o-prefix $case_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m FE >> {log.stdout} 2>> {log.stderr}
			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}_FE.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph
			sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph | bgzip -c > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.fold_enrichment_signal.bedgraph.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#

			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step2: BUILDING POISSON PVALUE SIGNAL TRACKS"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tAPPLICATION2: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tAPPLICATION3: ucsc-bedclip version:%s\\n" "366"
			printf "\\tAPPLICATION4: ucsc-bedgraphtobigwig version:%s\\n" "366"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_treat_pileup.bdg"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph"
			printf "\\tOUTPUT2: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph.gz"

			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s %s\\n" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_control_lambda.bdg --o-prefix $case_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m ppois"
			printf "\\t%s | %s\\n" "slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0" "bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph"
			printf "\\t%s\\n" "sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph | bgzip -c > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph.gz"
			printf "\\t%s\\n" "tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph.gz"
			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_Name}}_control_lambda.bdg --o-prefix $case_Name \
			--outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois >> {log.stdout} 2>> {log.stderr}
			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph
			sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph | bgzip -c > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_Name}}.poisson_pvalue_signal.bedgraph.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "#############################################################################"
			case_vs_control_Name=$(basename {input.narrowPeak})
			control_Name=$(echo "${{case_vs_control_Name}}" | sed -e "s/^${{case_Name}}_VS_//" -e "s/.narrowPeak.gz$//")
			printf "%s\\n" "#############################################################################"
			printf "%s\\n" "SIGNAL PROCESSING: ${{control_Name}}"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step1: BUILDING FOLD ENRICHMENT SIGNAL TRACKS"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tAPPLICATION2: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tAPPLICATION3: ucsc-bedclip version:%s\\n" "366"
			printf "\\tAPPLICATION4: ucsc-bedgraphtobigwig version:%s\\n" "366"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_treat_pileup.bdg"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph"
			printf "\\tOUTPUT2: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph.gz"

			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_control_lambda.bdg --o-prefix $control_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m FE "
			printf "\\t%s | %s\\n" "slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}_FE.bdg -g {REFDIR}/mm10.chrom.sizes -b 0" "bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph"
			printf "\\t%s\\n" "sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph | bgzip -c >  {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph.gz"
			printf "\\t%s\\n" "tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph.gz"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_control_lambda.bdg \
			--o-prefix $control_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m FE >> {log.stdout} 2>> {log.stderr}
			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}_FE.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph
			sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph | bgzip -c > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.fold_enrichment_signal.bedgraph.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#

			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step2: BUILDING POISSON PVALUE SIGNAL TRACKS"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tAPPLICATION2: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tAPPLICATION3: ucsc-bedclip version:%s\\n" "366"
			printf "\\tAPPLICATION4: ucsc-bedgraphtobigwig version:%s\\n" "366"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_treat_pileup.bdg"
			printf "\\tOUTPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph"
			printf "\\tOUTPUT2: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph.gz"

			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s %s\\n" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_control_lambda.bdg --o-prefix $control_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m ppois"
			printf "\\t%s | %s\\n" "slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0" "bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph"
			printf "\\t%s\\n" "sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph | bgzip -c > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph.gz"
			printf "\\t%s\\n" "tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph.gz"
			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{control_Name}}_control_lambda.bdg --o-prefix $control_Name \
			--outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois >> {log.stdout} 2>> {log.stderr}
			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph
			sort -k1,1 -k2,2n {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph | bgzip -c > {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph.gz
			tabix -f -p bed {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{control_Name}}.poisson_pvalue_signal.bedgraph.gz
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "#############################################################################"
			case_vs_control_Name=$(basename {input.narrowPeak})
			case_vs_control_Name=${{case_vs_control_Name/%.narrowPeak.gz/}}
			printf "%s\\n" "#############################################################################"
			printf "%s\\n" "SIGNAL PROCESSING: ${{case_Name}}_VS_${{control_Name}}"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step1: BUILDING FOLD ENRICHMENT SIGNAL TRACKS"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tAPPLICATION2: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tAPPLICATION3: ucsc-bedclip version:%s\\n" "366"
			printf "\\tAPPLICATION4: ucsc-bedgraphtobigwig version:%s\\n" "366"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_treat_pileup.bdg"
			printf "\\tOUTPUT1: %s\\n" "{output.fold_enrichment_signal_Bedgraph}"
			printf "\\tOUTPUT2: %s\\n" "{output.fold_enrichment_sorted_signal_Bedgraph_gz}"

			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s\\n" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_control_lambda.bdg --o-prefix $case_vs_control_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m FE "
			printf "\\t%s | %s\\n" "slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_vs_control_Name}}_FE.bdg -g {REFDIR}/mm10.chrom.sizes -b 0" "bedClip stdin {REFDIR}/mm10.chrom.sizes {output.fold_enrichment_signal_Bedgraph}"
			printf "\\t%s\\n" "sort -k1,1 -k2,2n {output.fold_enrichment_signal_Bedgraph} | bgzip -c > {output.fold_enrichment_sorted_signal_Bedgraph_gz}"
			printf "\\t%s\\n" "tabix -f -p bed {output.fold_enrichment_sorted_signal_Bedgraph_gz}"

			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_control_lambda.bdg \
			--o-prefix $case_vs_control_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m FE >> {log.stdout} 2>> {log.stderr}
			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_vs_control_Name}}_FE.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {output.fold_enrichment_signal_Bedgraph}
			sort -k1,1 -k2,2n {output.fold_enrichment_signal_Bedgraph} | bgzip -c > {output.fold_enrichment_sorted_signal_Bedgraph_gz}
			tabix -f -p bed {output.fold_enrichment_sorted_signal_Bedgraph_gz}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#

			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++"
			printf "%s\\n" "Step2: BUILDING POISSON PVALUE SIGNAL TRACKS"
			printf "\\tAPPLICATION1: macs version:%s\\n" "$(macs2 --version 2>&1 | head -n 1)"
			printf "\\tAPPLICATION2: bedtools version:%s\\n" "$(bedtools --version | head -1 | cut -f2 -d'v')"
			printf "\\tAPPLICATION3: ucsc-bedclip version:%s\\n" "366"
			printf "\\tAPPLICATION4: ucsc-bedgraphtobigwig version:%s\\n" "366"
			printf "\\tINPUT1: %s\\n" "{WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_treat_pileup.bdg"
			printf "\\tOUTPUT1: %s\\n" "{output.poisson_pvalue_signal_Bedgraph}"
			printf "\\tOUTPUT2: %s\\n" "{output.poisson_pvalue_sorted_signal_Bedgraph_gz}"

			printf "\\tCOMMAND_LINE:\\n"
			printf "\\t%s %s\\n" "macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_control_lambda.bdg --o-prefix $case_vs_control_Name --outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2 -m ppois"
			printf "\\t%s | %s\\n" "slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_vs_control_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0" "bedClip stdin {REFDIR}/mm10.chrom.sizes {output.poisson_pvalue_signal_Bedgraph}"
			printf "\\t%s\\n" "sort -k1,1 -k2,2n {output.poisson_pvalue_signal_Bedgraph} | bgzip -c > {output.poisson_pvalue_sorted_signal_Bedgraph_gz}"
			printf "\\t%s\\n" "tabix -p bed {output.poisson_pvalue_sorted_signal_Bedgraph_gz}"
			#
			printf "\\n"
			printf "\\t%s\\n" "EXECUTING...."
			start_time="$(date -u +%s)"
			macs2 bdgcmp -t {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_treat_pileup.bdg -c {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/${{case_vs_control_Name}}_control_lambda.bdg --o-prefix $case_vs_control_Name \
			--outdir {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/ -m ppois >> {log.stdout} 2>> {log.stderr}
			slopBed -i {WORKDIR}/{TITLE}/{wildcards.study}/SIGNAL_PROCESSING/MACS2/${{case_vs_control_Name}}_ppois.bdg -g {REFDIR}/mm10.chrom.sizes -b 0 | bedClip stdin {REFDIR}/mm10.chrom.sizes {output.poisson_pvalue_signal_Bedgraph}
			sort -k1,1 -k2,2n {output.poisson_pvalue_signal_Bedgraph} | bgzip -c > {output.poisson_pvalue_sorted_signal_Bedgraph_gz}
			tabix -p bed {output.poisson_pvalue_sorted_signal_Bedgraph_gz}
			end_time="$(date -u +%s)"
			printf "\\t%s\\n" "DONE!!!!"
			printf "\\tELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))"
			printf "\\n"
			#
			printf "%s\\n" "----------------------------------------"
			printf "%s\\n" "#############################################################################"
			""")
'''

