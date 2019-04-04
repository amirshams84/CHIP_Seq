shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for PreProcess
# snakemake --snakefile PreProcess.py --configfile Yoko.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile PreProcess.py --configfile Yoko.json --rulegraph | dot -Tsvg > CHIP_Seq.svg
# ################################### IMPORT ##################################### #


import os
import sys
import re
import glob
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


def get_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + SAMPLE_DELIMITER + "*" + "." + SAMPLE_SUFFIX, recursive=True)


def get_forward_fastq(wildcards):
	"""
	"""
	return glob.glob(DATADIR + "/**/" + wildcards.sample + SAMPLE_DELIMITER + "*" + FORWARD_DELIMITER + "*" + "." + SAMPLE_SUFFIX, recursive=True)
# ################################### CONFIGURATION ############################## #

# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
INFOLINK = config_general_Dict["INFOLINK"]
TITLE = config_general_Dict["TITLE"]
EXECUTION_MODE = config_general_Dict["EXECUTION_MODE"]
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
DATADIR = utility.fix_path(config_general_Dict["DATADIR"])
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
PLATFORM = config_data_Dict["PLATFORM"].lower()
FORMAT = config_data_Dict["FORMAT"].lower()
LAYOUT = config_data_Dict["LAYOUT"].lower()
SAMPLE_DELIMITER = config_data_Dict["SAMPLE_DELIMITER"]
SAMPLE_SUFFIX = config_data_Dict["SAMPLE_SUFFIX"].lower()
GENOME = config_data_Dict["GENOME"].lower()
####
if LAYOUT == "paired":
	#
	if PLATFORM == "illumina":
		FORWARD_DELIMITER = "R1"
		REVERSE_DELIMITER = "R2"
	elif PLATFORM == "sra":
		FORWARD_DELIMITER = ".1"
		REVERSE_DELIMITER = ".2"
	else:
		print("This pipeline only can be used with illumina-sra paired-end fastq format")
		print("Aborting!!!!")
		sys.exit(2)
elif LAYOUT == "single":
	FORWARD_DELIMITER = ""
	REVERSE_DELIMITER = ""
else:
	print("This pipeline only can be used with single/paired-end fastq format")
	print("Aborting!!!!")
	sys.exit(2)
####
if SAMPLE_SUFFIX[0] == ".":
	SAMPLE_SUFFIX = SAMPLE_SUFFIX[1:]
if SAMPLE_SUFFIX not in ["fastq", "fastq.gz"]:
	print("This pipeline only can bed used with fastq or fastq.gz data")
	print("Aborting!!!!")
	sys.exit(2)
else:
	pass
####
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
design_Dict = build_design_Dict(metadata_Dict)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#PRE_PROCESS
config_pre_process_Dict = config["PRE_PROCESS"]
# ------------------------------------
# ################################### WILDCARDS ################################ #


pre_process_List = []
for sample, sample_Dict in metadata_Dict.items():
	pre_process_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq.gz".format(design=sample_Dict["Design"], sample=sample))
# ################################### PIPELINE FLOW ############################ #

rule End_Point:
	input:
		pre_process_List
# ################################### PIPELINE RULES ########################## #

if LAYOUT == "paired":
	#
	rule PreProcess_Paired:
		"""
		"""
		input:
			fastq_List = get_forward_fastq
		output:
			processed_fwd_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq.gz",
			processed_rev_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R2.processed.fastq.gz",
		threads: PROCESSORS
		message: "PreProcess_Paired: {wildcards.design}|{wildcards.sample}"
		resources:
			mem_mb = MEMORY
		run:
			for each_fastq in input.fastq_List:
				#
				each_fastq_reverse = each_fastq.replace(FORWARD_DELIMITER, REVERSE_DELIMITER)
				each_fastq_basename = os.path.basename(each_fastq)
				each_fastq_begining = re.sub("." + SAMPLE_SUFFIX, "", each_fastq_basename)
				shell("""
					#
					module load cutadapt/2.1
					QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/pre_process
					mkdir -p $QC_PATH
					#
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "CUTADAPT2.1"  | tee >(cat >&2)
					printf "%s\\n" "Removing low complexity, remnant adapter and short reads."  | tee >(cat >&2)
					printf "INPUT1: %s\\n" "{each_fastq}"  | tee >(cat >&2)
					printf "INPUT2: %s\\n" "{each_fastq_reverse}"  | tee >(cat >&2)
					printf "OUTPUT1: %s\\n" "{output.processed_fwd_fastq}.tmp"  | tee >(cat >&2)
					printf "OUTPUT2: %s\\n" "{output.processed_rev_fastq}.tmp"  | tee >(cat >&2)
					printf "OUTPUT3: %s\\n" "$QC_PATH/{each_fastq_begining}.cutadapt.txt"  | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
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
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "CAT"  | tee >(cat >&2)
					printf "%s\\n" "Aggregating processed reads from technical replicates into single sample"  | tee >(cat >&2)
					printf "INPUT1: %s\\n" "{output.processed_fwd_fastq}.tmp"  | tee >(cat >&2)
					printf "INPUT2: %s\\n" "{output.processed_rev_fastq}.tmp"  | tee >(cat >&2)
					printf "OUTPUT1: %s\\n" "{output.processed_fwd_fastq}"  | tee >(cat >&2)
					printf "OUTPUT2: %s\\n" "{output.processed_rev_fastq}"  | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
					printf "%s\\n" "cat {output.processed_fwd_fastq}.tmp >> {output.processed_fwd_fastq}" | tee >(cat >&2)
					printf "%s\\n" "cat {output.processed_rev_fastq}.tmp >> {output.processed_rev_fastq}" | tee >(cat >&2)
					printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
					start_time="$(date -u +%s)"
					#
					##
					cat {output.processed_fwd_fastq}.tmp >> {output.processed_fwd_fastq}
					cat {output.processed_rev_fastq}.tmp >> {output.processed_rev_fastq}
					rm {output.processed_fwd_fastq}.tmp {output.processed_rev_fastq}.tmp
					##
					#
					end_time="$(date -u +%s)"
					printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
				""")
elif LAYOUT == "single":
	#
	rule PreProcess_Single:
		input:
			fastq_List = get_fastq
		output:
			processed_fastq = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/pre_process/{sample}.R1.processed.fastq.gz",
		threads: PROCESSORS
		message: "PreProcess_Single: {wildcards.design}|{wildcards.sample}"
		resources:
			mem_mb = MEMORY
		run:
			for each_fastq in input.fastq_List:
				#
				each_fastq_basename = os.path.basename(each_fastq)
				each_fastq_begining = re.sub("." + SAMPLE_SUFFIX, "", each_fastq_basename)
				shell("""
					module load cutadapt/2.1
					QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/pre_process
					mkdir -p $QC_PATH
					#
					printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
					printf "%s\\n" "CUTADAPT2.1"  | tee >(cat >&2)
					printf "%s\\n" "Removing low complexity, remnant adapter and short reads."  | tee >(cat >&2)
					printf "INPUT1: %s\\n" "{each_fastq}"  | tee >(cat >&2)
					printf "OUTPUT1: %s\\n" "{output.processed_fastq}"  | tee >(cat >&2)
					printf "OUTPUT3: %s\\n" "$QC_PATH/{each_fastq_begining}.cutadapt.txt"  | tee >(cat >&2)
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
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
					printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)

				""")
