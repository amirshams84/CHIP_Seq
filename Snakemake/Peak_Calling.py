shell.executable("/bin/bash")
shell.prefix("source /data/shamsaddinisha/conda/etc/profile.d/conda.sh")
# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Peak Calling
# snakemake --snakefile Peak_Calling.py --configfile Yoko.json --cores=50 -j 10 --local-cores=10
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


post_alignment_List = []
peak_calling_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{sample}.processed.samtools.txt".format(design=sample_Dict["Design"], sample=sample))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))


for design in design_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/report/post_alignment/{pooled_case}.processed.samtools.txt".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##
	for case in design_Dict[design]["Case"]:
		for control in design_Dict[design]["Control"]:
			#
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
	##
	for control in design_Dict[design]["Control"]:
		#
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{pooled_case}_VS_{control}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		peak_calling_List
# ################################### PIPELINE RULES ########################## #

#+++++++++++++++++++++++++++++
##REALM4: PEAK-CALLING
#+++++++++++++++++++++++++++++

rule Peak_Calling_Narrow:
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_index_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
	output:
		#
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.narrowPeak.gz",
		narrowPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.narrowPeak.bb",
		narrowPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.narrowPeak.FE.bigwig",
		narrowPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.narrowPeak.PV.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Narrow: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {input.processed_bam})
			sample_Name=${{sample_Name%.processed.bam}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling
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
			if [ -f {output.narrowPeak_PV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg
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
		processed_case_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
	output:
		#
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.gz",
		narrowPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.bb",
		narrowPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.FE.bigwig",
		narrowPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.narrowPeak.PV.bigwig",
		
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Narrow_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	run:
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {output.narrowPeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling
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
			if [ -f {output.narrowPeak_PV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_narrow_ppois.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.FE.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.narrowPeak.PV.bigwig.bdg
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
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_index_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
	output:
		#
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.broadPeak.gz",
		broadPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.broadPeak.bb",
		broadPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.broadPeak.FE.bigwig",
		broadPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{sample}.broadPeak.PV.bigwig",
		
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {input.processed_bam})
			sample_Name=${{sample_Name%.processed.bam}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling
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
			if [ -f {output.broadPeak_PV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
				start_time="$(date -u +%s)"
				#
				##
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_qpois.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg
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
		processed_case_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_index_Bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
	output:
		#
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.gz",
		broadPeak_index_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.bb",
		broadPeak_FE_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.FE.bigwig",
		broadPeak_PV_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/{case}_VS_{control}.broadPeak.PV.bigwig",
		
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	run:
		shell("""
			echo "{ACTIVATE_CONDA_PY2}"
			{ACTIVATE_CONDA_PY2}
			sample_Name=$(basename {output.broadPeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling
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
			if [ -f {output.broadPeak_PV_bigwig} ]; then

				printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg $OUT_PATH/${{sample_Name}}.macs2_broad_ppois.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg" | tee >(cat >&2)
				printf "%s\\n" "rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg" | tee >(cat >&2)
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
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.FE.bigwig.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg.tmp $OUT_PATH/${{sample_Name}}.broadPeak.PV.bigwig.bdg

				##
				#
				end_time="$(date -u +%s)"
				printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			fi
		""")