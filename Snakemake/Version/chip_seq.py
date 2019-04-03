shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Aug-14-2018
# Email: amir.shams84@gmail.com
# Aim: Snakemake recepie for CHIP-Seq Comprehensive analysis
# snakemake --snakefile chip_seq.py --configfile chip_seq.json --debug-dag --cores=50
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


def get_case_bam(wildcards):
	processed_bam_List = []
	for Index, study_Dict in config_metadata_Dict.items():
		if study_Dict["Title"] == wildcards.study:
			for case in study_Dict["Case"]:
				processed_bam_List.append(WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{case}.processed.bam".format(study=wildcards.study, case=case))
	return processed_bam_List


def get_narrowPeak(wildcards):
	narrowPeak_List = []
	for Index, study_Dict in config_metadata_Dict.items():
		if study_Dict["Title"] == wildcards.study:
			for case in study_Dict["Case"]:
				narrowPeak_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{case}.narrowPeak.gz".format(study=wildcards.study, case=case))
	return narrowPeak_List


def get_narrowPeak_pooled(wildcards):
	narrowPeak_List = []
	for Index, study_Dict in config_metadata_Dict.items():
		if study_Dict["Title"] == wildcards.study:
			for case in study_Dict["Case"]:
				narrowPeak_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{case}_VS_{control}.narrowPeak.gz".format(study=study_Dict["Title"], case=case, control=control))
	
	return narrowPeak_List


# ################################### CONFIGURATION ############################## #
#((?!_VS_).)*

localrules: all
configfile: "chip_seq_pooling_idr.json"
config_general_Dict = config["GENERAL"]
config_data_Dict = config["DATA"]
config_metadata_Dict = config["METADATA"]
config_pre_processing_Dict = config["PRE_PROCESSING"]
config_alignment_Dict = config["ALIGNMENT"]
config_post_alignment_Dict = config["POST_ALIGNMENT"]
config_peak_calling_Dict = config["PEAK_CALLING"]
config_reference_Dict = config["REFERENCE"]
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

POOLED_CASE_List = []
POOLED_NARROW_PEAK_CALLING_List = []

IDR_List = []
for Index, study_Dict in config_metadata_Dict.items():
	#
	IDR_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case}.IDR.txt".format(study=study_Dict["Title"], pooled_case="_POOLED_".join(study_Dict["Case"])))
	POOLED_CASE_List.append(WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{pooled_case}.processed.bam".format(study=study_Dict["Title"], pooled_case="_POOLED_".join(study_Dict["Case"])))
	POOLED_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case}.narrowPeak.gz".format(study=study_Dict["Title"], pooled_case="_POOLED_".join(study_Dict["Case"])))
	for case in study_Dict["Case"]:
		#
		if case not in case_List:
			case_List.append("{case}".format(case=case))
			case_PRE_PROCESSING_List.append(WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{case}.R1.processed.fastq.gz".format(study=study_Dict["Title"], case=case))
			case_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{case}.bam".format(study=study_Dict["Title"], case=case))
			case_POST_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{case}.processed.bam".format(study=study_Dict["Title"], case=case))
			case_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{case}.narrowPeak.gz".format(study=study_Dict["Title"], case=case))
	for control in study_Dict["Control"]:
		#
		if control not in control_List:
			control_List.append("{control}".format(control=control))
			control_PRE_PROCESSING_List.append(WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{control}.R1.processed.fastq.gz".format(study=study_Dict["Title"], control=control))
			control_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{control}.bam".format(study=study_Dict["Title"], control=control))
			control_POST_ALIGNMENT_List.append(WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{control}.processed.bam".format(study=study_Dict["Title"], control=control))
			control_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{control}.narrowPeak.gz".format(study=study_Dict["Title"], control=control))


sample_List = case_List + control_List
PRE_PROCESSING_List = case_PRE_PROCESSING_List + control_PRE_PROCESSING_List
ALIGNMENT_List = case_ALIGNMENT_List + control_ALIGNMENT_List
POST_ALIGNMENT_List = case_POST_ALIGNMENT_List + control_POST_ALIGNMENT_List
NARROW_PEAK_CALLING_List = case_NARROW_PEAK_CALLING_List + control_NARROW_PEAK_CALLING_List


controlled_NARROW_PEAK_CALLING_List = []
controlled_POOLED_NARROW_PEAK_CALLING_List = []
controlled_IDR_List = []
for Index, study_Dict in config_metadata_Dict.items():
	#
	for case in study_Dict["Case"]:
		#
		for control in study_Dict["Control"]:
			#
			controlled_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{case}_VS_{control}.narrowPeak.gz".format(study=study_Dict["Title"], case=case, control=control))
			controlled_POOLED_NARROW_PEAK_CALLING_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case}_VS_{control}.narrowPeak.gz".format(study=study_Dict["Title"], pooled_case="_POOLED_".join(study_Dict["Case"]), control=control))
			controlled_IDR_List.append(WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case}_VS_{control}.IDR.txt".format(study=study_Dict["Title"], pooled_case="_POOLED_".join(study_Dict["Case"]), control=control))
# ################################### RULES ###################################### #
rule all:
	input: PRE_PROCESSING_List + ALIGNMENT_List + POST_ALIGNMENT_List + NARROW_PEAK_CALLING_List + controlled_NARROW_PEAK_CALLING_List + POOLED_CASE_List + POOLED_NARROW_PEAK_CALLING_List + controlled_POOLED_NARROW_PEAK_CALLING_List + IDR_List + controlled_IDR_List

rule PRE_PROCESSING:
	input:
		fq1_List = get_fq1,
	output:
		processed_R1 = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample}.R2.processed.fastq.gz",
		discarded_R1 = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample}.R1.discarded.fastq.gz",
		discarded_R2 = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample}.R2.discarded.fastq.gz",
		cutadapt_report = WORKDIR + "/" + TITLE + "/{study}/QUALITY_CONTROL/{sample}.cutadapt_report.txt",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample}.PRE_PROCESSING.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample}.PRE_PROCESSING.stderr",
	priority: 99
	threads: PROCESSORS
	message: "PRE_PROCESSING"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			declare -a fq1_List=({input.fq1_List})
			for i in "${{fq1_List[@]}}"
			do
				fq1=$i
				fq2=${{fq1/%_R1_001.fastq.gz/_R2_001.fastq.gz}}
				fq1_Name=$(basename $fq1)
				trim_fq1={WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/${{fq1_Name/%_R1_001.fastq.gz/.R1.processed.fastq.gz}}
				trim_fq2={WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/${{fq1_Name/%_R1_001.fastq.gz/.R2.processed.fastq.gz}}
				discarded_fq1={WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/${{fq1_Name/%_R1_001.fastq.gz/.R1.discarded.fastq.gz}}
				discarded_fq2={WORKDIR}/{TITLE}/{wildcards.study}/PRE_PROCESSING/${{fq1_Name/%_R1_001.fastq.gz/.R2.discarded.fastq.gz}}
				echo "cutadapt {config_pre_processing_Dict[CUTADAPT]} --too-short-output=$discarded_fq1 --too-short-paired-output=$discarded_fq2 --output=$trim_fq1 --paired-output=$trim_fq2 $fq1 $fq2"
				cutadapt {config_pre_processing_Dict[CUTADAPT]} --too-short-output=$discarded_fq1 --too-short-paired-output=$discarded_fq2 --output=$trim_fq1 --paired-output=$trim_fq2 $fq1 $fq2 > {output.cutadapt_report}
				zcat $discarded_fq1 | gzip -nc >> {output.discarded_R1}
				zcat $discarded_fq2 | gzip -nc >> {output.discarded_R2}
				zcat $trim_fq1 | gzip -nc >> {output.processed_R1}
				zcat $trim_fq2 | gzip -nc >> {output.processed_R2}
			done
			""")


rule ALIGNMENT:
	input:
		processed_R1 = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample, [^/.]+}.R1.processed.fastq.gz",
		processed_R2 = WORKDIR + "/" + TITLE + "/{study}/PRE_PROCESSING/{sample, [^/.]+}.R2.processed.fastq.gz"
	output:
		alignment_Bam = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, [^/.]+}.bam",
		alignment_Report = WORKDIR + "/" + TITLE + "/{study}/QUALITY_CONTROL/{sample, [^/.]+}.alignment_report.txt",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, [^/.]+}.ALIGNMENT.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, [^/.]+}.ALIGNMENT.stderr",

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
		alignment_Bam = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.bam",
	output:
		#PICARD
		marked_duplicate_Bam = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.marked_duplicate.bam",
		picard_report = WORKDIR + "/" + TITLE + "/{study}/QUALITY_CONTROL/{sample, ((?!_POOLED_).)*}.picard_report.txt",
		flagstat = WORKDIR + "/" + TITLE + "/{study}/QUALITY_CONTROL/{sample, ((?!_POOLED_).)*}.flagstat.txt",
		#SAMTOOLS
		processed_Bam = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.processed.bam",
		processed_Flagstat = WORKDIR + "/" + TITLE + "/{study}/QUALITY_CONTROL/{sample, ((?!_POOLED_).)*}.processed.flagstat.txt",
		nmsort_Bam = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.processed.nmsort.bam",
		Bed = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.bed.gz",
		processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.processed.tagAlign.gz",
	log:
		stdout = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.POST_ALIGNMENT.stdout",
		stderr = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_POOLED_).)*}.POST_ALIGNMENT.stderr",
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


rule NARROW_PEAK_CALLING:
	input:
		processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{sample, ((?!_VS_).)*}.processed.bam",
	output:
		narrowPeak = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{sample, ((?!_VS_).)*}.narrowPeak.gz",
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
			sample_Name=${{sample_Name%.processed.bam}}
			echo "macs2 callpeak --treatment {input.processed_tagAlign} -f BAMPE --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			macs2 callpeak --treatment {input.processed_tagAlign} -f BAMPE --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/
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


rule NARROW_PEAK_CALLING_CONTROLLED:
	input:
		case_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{case}.processed.bam",
		control_processed_tagAlign = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{control}.processed.bam",
	output:
		versus_narrowPeak = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{case}_VS_{control}.narrowPeak.gz",
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
			echo "macs2 callpeak --treatment {input.case_processed_tagAlign} --control {input.control_processed_tagAlign}  -f BAMPE --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/"
			macs2 callpeak --treatment {input.case_processed_tagAlign} --control {input.control_processed_tagAlign}  -f BAMPE --name $sample_Name {config_peak_calling_Dict[MACS2_NARROW]} --outdir {WORKDIR}/{TITLE}/{wildcards.study}/PEAK_CALLING/MACS2/
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


rule POOLING_CASE_REPLICATE:
	input:
		case_bam_List = get_case_bam
	output:
		pooled_bam = WORKDIR + "/" + TITLE + "/{study}/ALIGNMENT/{pooled_case, .*_POOLED_.*}.processed.bam",
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY2}
			echo "pooling cases"
			echo "Pooling: {input.case_bam_List}"
			echo "samtools merge {output.pooled_bam} {input.case_bam_List}"
			samtools merge {output.pooled_bam} {input.case_bam_List}
			samtools index {output.pooled_bam}
			""")


rule IDR:
	input:
		narrowPeak_List = get_narrowPeak,
		pooled_narrowPeak = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case, ((?!_VS_).)*}.narrowPeak.gz",
	output:
		idr_result = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case, ((?!_VS_).)*}.IDR.txt",
		idr_log = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case, ((?!_VS_).)*}.IDR.log",
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY3}
			printf "\\t%s\\n" "idr --samples {input.narrowPeak_List} --peak-list {input.pooled_narrowPeak} --output-file {output.idr_result}  --log-output-file {output.idr_log}"

			idr --samples {input.narrowPeak_List} --peak-list {input.pooled_narrowPeak} --output-file {output.idr_result}  --log-output-file {output.idr_log}
			""")

rule IDR_CONTROLLED:
	input:
		narrowPeak_List = get_narrowPeak_pooled,
		pooled_narrowPeak = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case}_VS_{control}.narrowPeak.gz",
	output:
		idr_result = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case}_VS_{control}.IDR.txt",
		idr_log = WORKDIR + "/" + TITLE + "/{study}/PEAK_CALLING/MACS2/{pooled_case}_VS_{control}.IDR.log",
	run:
		shell("""
			{CONDA_INIT}
			{ACTIVATE_ENV_PY3}
			printf "\\t%s\\n" "idr --samples {input.narrowPeak_List} --peak-list {input.pooled_narrowPeak} --output-file {output.idr_result}  --log-output-file {output.idr_log}"

			idr --samples {input.narrowPeak_List} --peak-list {input.pooled_narrowPeak} --output-file {output.idr_result}  --log-output-file {output.idr_log}
			""")