shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Alignment
# snakemake --snakefile PreProcess.py --configfile Encode.json --cores=50 -j 1 --local-cores=10
# snakemake --snakefile PreProcess.py --configfile Yoko.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile PreProcess.py --configfile CHIP_Seq.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
# snakemake --snakefile PreProcess.py --configfile Yoko.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
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
	pre_process_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam".format(design=sample_Dict["Design"], sample=sample))
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{sample}.samtools.txt".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	#
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case}.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{pooled_case}.samtools.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
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
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq.gz",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R2.processed.fastq.gz",
		output:
			bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
			bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
			bowtie2_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bigwig"
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "Alignment_Paired: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fwd_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq.gz", "", each_fastq_basename)
			shell("""
				echo "{ACTIVATE_CONDA_PY2}"
				{ACTIVATE_CONDA_PY2}
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "BOWTIE2:" | tee >(cat >&2)
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
				printf "%s\\n" "PICARD:" | tee >(cat >&2)
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
				printf "%s\\n" "SAMTOOLS:" | tee >(cat >&2)
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
				printf "%s\\n" "DEEPTOOLS:" | tee >(cat >&2)
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

elif LAYOUT == "single":
	rule Alignment_Single:
		"""
		"""
		input:
			processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq.gz",
		output:
			bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
			bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
			bowtie2_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bigwig",
		priority: 996
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		message: "Alignment_Single: {wildcards.design}|{wildcards.sample}"
		run:
			each_fastq_basename = os.path.basename(input.processed_fastq)
			each_fastq_begining = re.sub(".R1.processed.fastq.gz", "", each_fastq_basename)
			shell("""
				echo "{ACTIVATE_CONDA_PY2}"
				{ACTIVATE_CONDA_PY2}
				QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
				mkdir -p $QC_PATH
				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "BOWTIE2:" | tee >(cat >&2)
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
				printf "%s\\n" "PICARD:" | tee >(cat >&2)
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
				printf "%s\\n" "SAMTOOLS:" | tee >(cat >&2)
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
				printf "%s\\n" "DEEPTOOLS:" | tee >(cat >&2)
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
		pooled_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam",
		pooled_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bam.bai",
		pooled_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{pooled_case, .*_POOLED_.*}.bigwig"
	priority: 993
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Pooling_Case_Replicates: {wildcards.design}|{wildcards.pooled_case}"
	run:
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "SAMTOOLS:" | tee >(cat >&2)
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
			printf "%s\\n" "DEEPTOOLS:" | tee >(cat >&2)
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
		bowtie2_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam",
		bowtie2_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/alignment/{sample}.bam.bai",
	output:
		flagstats = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{sample}.samtools.txt",
		spp_phantomPeak = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{sample}.spp.pdf",
		spp_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{sample}.spp.txt",
		complexity_report = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/alignment/{sample}.complexity.txt",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "ALIGNMENT_QC: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/alignment
			mkdir -p $QC_PATH
			unset DISPLAY
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "QUALIMAP:" | tee >(cat >&2)
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
			printf "%s\\n" "FASTQC:" | tee >(cat >&2)
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
			printf "%s\\n" "SAMTOOLS:" | tee >(cat >&2)
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
			printf "%s\\n" "SPP:" | tee >(cat >&2)
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
