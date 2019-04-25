# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Mar-29-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for Peak Calling
# snakemake --snakefile peak_calling.py --configfile Encode.json --cores=50 -j 10 --local-cores=10
# snakemake --snakefile peak_calling.py --configfile Encode.json --rulegraph | dot -Tsvg > peak_calling.svg
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
MACS2_NARROW_PARAMETERS = config_peak_calling_Dict["MACS2_NARROW"]
MACS2_BROAD_PARAMETERS = config_peak_calling_Dict["MACS2_BROAD"]

# ------------------------------------
# ################################### WILDCARDS ################################ #


post_alignment_List = []
peak_calling_List = []
overlap_peak_List = []
for sample, sample_Dict in metadata_Dict.items():
	#
	#POST_ALIGNMENT
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam".format(design=sample_Dict["Design"], sample=sample))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz".format(design=sample_Dict["Design"], sample=sample))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz".format(design=sample_Dict["Design"], sample=sample))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz".format(design=sample_Dict["Design"], sample=sample))

for design in design_Dict:
	##POOLING
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bam".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	post_alignment_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{pooled_case}.processed.bed.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##PEAK_CALLING
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"])))
	##
	for case in design_Dict[design]["Case"]:
		for control in design_Dict[design]["Control"]:
			#
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz".format(design=design, case=case, control=control))
			peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz".format(design=design, case=case, control=control))
	##
	for control in design_Dict[design]["Control"]:
		#
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{pooled_case}_VS_{control}.narrowPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		peak_calling_List.append(WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{pooled_case}_VS_{control}.broadPeak.gz".format(design=design, pooled_case="_POOLED_".join(design_Dict[design]["Case"]), control=control))
		

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
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
		processed_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz",
		processed_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz.tbi",
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{sample}.narrowPeak.bigwig",
		
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Narrow: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			#
			module load samtools/1.9
			module load bedtools/2.27.1
			module load macs/2.1.2
			module load ucsc/373
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling
			mkdir -p $QC_PATH
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak
			mkdir -p $OUT_PATH
			#
			if [ ! -f ./Script/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O ./Script/bigNarrowPeak.as
			fi
			#
			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}
			#
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "narrow peak calling"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp
			bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "building signal track using Fold Enrichment"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg" | tee >(cat >&2)
			printf "INPUT2: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}
			#
			cut -f1,2,3 {output.narrowPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,2 {output.narrowPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,3 {output.narrowPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			if [ -s $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then
				while IFS= read line; do grep -nr -m 1 --perl-regex "$line" {output.narrowPeak_bdg} | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
				while IFS= read line; do sed -i "${{line}}d" {output.narrowPeak_bdg}; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
			fi
			#
			bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			if [ -f {output.narrowPeak_bigwig} ]; then
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bdg.tmp
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
			fi

		""")

rule Peak_Calling_Narrow_Controlled:
	input:
		processed_case_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
		#
		processed_case_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz",
		processed_case_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz.tbi",
		processed_control_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz",
		processed_control_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz.tbi",
	output:
		narrowPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz",
		narrowPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.gz.tbi",
		narrowPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bb",
		narrowPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bdg",
		narrowPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/narrowpeak/{case}_VS_{control}.narrowPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Narrow_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	run:
		shell("""
			#
			module load samtools/1.9
			module load bedtools/2.27.1
			module load macs/2.1.2
			module load ucsc/373
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/narrowpeak/
			mkdir -p $OUT_PATH
			#
			if [ ! -f ./Script/bigNarrowPeak.as ]; then
				wget {config_utilities_Dict[BigNarrowPeak]} -O ./Script/bigNarrowPeak.as
			fi
			#
			sample_Name=$(basename {output.narrowPeak_bed})
			sample_Name=${{sample_Name%.narrowPeak.gz}}
			#
			AWK_COMMAND='BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "narrow peak calling"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_case_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_control_bam}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_case_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_control_bed}"  | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.narrowPeak_bigbed}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.narrowPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_narrow {MACS2_NARROW_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak > {output.narrowPeak_bed}.sorted
			bgzip -c {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}
			tabix -f -p bed {output.narrowPeak_bed}
			awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.narrowPeak_bed}.sorted > {output.narrowPeak_bed}.tmp
			bedToBigBed -as=./Script/bigNarrowPeak.as -type=bed6+4 {output.narrowPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigbed}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "building signal track using Fold Enrichment"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg" | tee >(cat >&2)
			printf "INPUT2: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.narrowPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_narrow --outdir $OUT_PATH --method FE
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.narrowPeak_bdg}.tmp > {output.narrowPeak_bdg}
			#
			cut -f1,2,3 {output.narrowPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,2 {output.narrowPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			cut -f1,3 {output.narrowPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
			if [ -s $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt ]; then
				while IFS= read line; do grep -nr -m 1 --perl-regex "$line" {output.narrowPeak_bdg} | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
				while IFS= read line; do sed -i "${{line}}d" {output.narrowPeak_bdg}; done< $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
			fi
			#
			bedGraphToBigWig {output.narrowPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.narrowPeak_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			if [ -f {output.narrowPeak_bigwig} ]; then
				#rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.xls
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_peaks.narrowPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_narrow_FE.bdg
				rm -rf {output.narrowPeak_bed}.tmp
				rm -rf {output.narrowPeak_bed}.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.narrowPeak.bdg.tmp
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.narrow.duplicate_line.txt
			fi
		""")


rule Peak_Calling_Broad:
	input:
		processed_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam",
		processed_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bam.bai",
		processed_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz",
		processed_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{sample}.processed.bed.gz.tbi",
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{sample}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad: {wildcards.design}|{wildcards.sample}"
	run:
		shell("""
			#
			module load samtools/1.9
			module load bedtools/2.27.1
			module load macs/2.1.2
			module load ucsc/373
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling
			mkdir -p $QC_PATH
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $OUT_PATH
			#
			if [ ! -f ./Script/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O ./Script/bigBroadPeak.as
			fi
			#
			sample_Name=$(basename {input.processed_bed})
			sample_Name=${{sample_Name%.processed.bed.gz}}
			#
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "Broad peak calling"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_bam_index}"  | tee >(cat >&2)
			printf "INPUT3: %s\\n" "{input.processed_bed}"  | tee >(cat >&2)
			printf "INPUT4: %s\\n" "{input.processed_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp
			bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "building signal track"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg" | tee >(cat >&2)
			printf "INPUT2: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}
			#
			cut -f1,2,3 {output.broadPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,2 {output.broadPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,3 {output.broadPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			if [ -s $OUT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then
				while IFS= read line; do grep -nr -m 1 --perl-regex "$line" {output.broadPeak_bdg} | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
				while IFS= read line; do sed -i "${{line}}d" {output.broadPeak_bdg}; done< $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt
			fi
			#
			bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			if [ -f {output.broadPeak_bigwig} ]; then
				#rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bdg.tmp
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt

			fi
		""")

rule Peak_Calling_Broad_Controlled:
	input:
		processed_case_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam",
		processed_case_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bam.bai",
		processed_control_bam = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam",
		processed_control_bam_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bam.bai",
		#
		processed_case_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz",
		processed_case_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{case}.processed.bed.gz.tbi",
		processed_control_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz",
		processed_control_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/post_alignment/{control}.processed.bed.gz.tbi",
	output:
		broadPeak_bed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz",
		broadPeak_bed_index = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.gz.tbi",
		broadPeak_bigbed = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bb",
		broadPeak_bdg = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bdg",
		broadPeak_bigwig = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE + "/" + GENOME + "/{design}/peak_calling/broadpeak/{case}_VS_{control}.broadPeak.bigwig",
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	message: "Peak_Calling_Broad_Controlled: {wildcards.design}|{wildcards.case}_VS_{wildcards.control}"
	run:
		shell("""
			#
			module load samtools/1.9
			module load bedtools/2.27.1
			module load macs/2.1.2
			module load ucsc/373
			QC_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/report/peak_calling
			mkdir -p $QC_PATH
			OUT_PATH={WORKDIR}/{PROJECT}/{EXPERIMENT}/{TITLE}/{GENOME}/{wildcards.design}/peak_calling/broadpeak
			mkdir -p $OUT_PATH
			#
			if [ ! -f ./Script/bigBroadPeak.as ]; then
				wget {config_utilities_Dict[BigBroadPeak]} -O ./Script/bigBroadPeak.as
			fi
			#
			sample_Name=$(basename {output.broadPeak_bed})
			sample_Name=${{sample_Name%.broadPeak.gz}}
			#
			AWK_COMMAND='BEGIN{{OFS="\\t"}}{{if ($5>1000) $5=1000; print $0}}'
			#
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "Broad peak calling"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_case_bam}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_control_bam}"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "{input.processed_case_bed}"  | tee >(cat >&2)
			printf "INPUT2: %s\\n" "{input.processed_control_bed}"  | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bed}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bed_index}"  | tee >(cat >&2)
			printf "OUTPUT3: %s\\n" "{output.broadPeak_bigbed}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH" | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.narrowPeak > {output.broadPeak_bed}.sorted" | tee >(cat >&2)
			printf "%s\\n" "bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}"  | tee >(cat >&2)
			printf "%s\\n" "tabix -f -p bed {output.broadPeak_bed}" | tee >(cat >&2)
			printf "%s\\n" "awk '$AWK_COMMAND' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp" | tee >(cat >&2)
			printf "%s\\n" "bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}" | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 callpeak --treatment {input.processed_case_bed} --control {input.processed_control_bed} --name ${{sample_Name}}.macs2_broad {MACS2_BROAD_PARAMETERS} --outdir $OUT_PATH
			LC_COLLATE=C sort -k1,1 -k2,2n $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak > {output.broadPeak_bed}.sorted
			bgzip -c {output.broadPeak_bed}.sorted > {output.broadPeak_bed}
			tabix -f -p bed {output.broadPeak_bed}
			awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' {output.broadPeak_bed}.sorted > {output.broadPeak_bed}.tmp
			bedToBigBed -as=./Script/bigBroadPeak.as -type=bed6+3 {output.broadPeak_bed}.tmp {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigbed}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			printf "%s\\n" "###################################- COMMANDLINE -############################" | tee >(cat >&2)
			printf "%s\\n" "samtools/1.9" | tee >(cat >&2)
			printf "%s\\n" "bedtools/2.27.1" | tee >(cat >&2)
			printf "%s\\n" "macs/2.1.2" | tee >(cat >&2)
			printf "%s\\n" "ucsc/373" | tee >(cat >&2)
			printf "%s\\n" "building signal track"  | tee >(cat >&2)
			printf "INPUT1: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg" | tee >(cat >&2)
			printf "INPUT2: %s\\n" "$OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg" | tee >(cat >&2)
			printf "OUTPUT1: %s\\n" "{output.broadPeak_bdg}"  | tee >(cat >&2)
			printf "OUTPUT2: %s\\n" "{output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" | tee >(cat >&2)
			printf "%s\\n" "macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE" | tee >(cat >&2)
			printf "%s\\n" "slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp"  | tee >(cat >&2)
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}" | tee >(cat >&2)
			printf "%s\\n" "bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}"  | tee >(cat >&2)
			printf "%s\\n" "EXECUTING...." | tee >(cat >&2)
			start_time="$(date -u +%s)"
			#
			##
			macs2 bdgcmp -t $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg -c $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg --o-prefix ${{sample_Name}}.macs2_broad --outdir $OUT_PATH --method FE
			slopBed -i $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg -g {config_reference_Dict[CHROM_SIZE]} -b 0 | bedClip stdin {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bdg}.tmp
			LC_COLLATE=C sort -k1,1 -k2,2n {output.broadPeak_bdg}.tmp > {output.broadPeak_bdg}
			#
			cut -f1,2,3 {output.broadPeak_bdg} | sort | uniq -d > $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,2 {output.broadPeak_bdg} | sort | uniq -d >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			cut -f1,3 {output.broadPeak_bdg} | sort | uniq -d | sed 's/\\t/\\t.*\\t/' >> $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
			if [ -s $OUT_PATH/${{sample_Name}}.broad.duplicate.txt ]; then
				while IFS= read line; do grep -nr -m 1 --perl-regex "$line" {output.broadPeak_bdg} | cut -d":" -f1 >> $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt; done< $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
				while IFS= read line; do sed -i "${{line}}d" {output.broadPeak_bdg}; done< $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt
			fi
			#
			bedGraphToBigWig {output.broadPeak_bdg} {config_reference_Dict[CHROM_SIZE]} {output.broadPeak_bigwig}
			##
			#
			end_time="$(date -u +%s)"
			printf "%s\\n" "DONE!!!!" | tee >(cat >&2)
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" | tee >(cat >&2)
			printf "%s\\n" "------------------------------------------------------------------------------" | tee >(cat >&2)
			if [ -f {output.broadPeak_bigwig} ]; then
				#rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.xls
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.broadPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_peaks.gappedPeak
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_summits.bed
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_treat_pileup.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_control_lambda.bdg
				rm -rf $OUT_PATH/${{sample_Name}}.macs2_broad_FE.bdg
				rm -rf {output.broadPeak_bed}.tmp
				rm -rf {output.broadPeak_bed}.sorted
				rm -rf $OUT_PATH/${{sample_Name}}.broadPeak.bdg.tmp
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate.txt
				rm -rf $OUT_PATH/${{sample_Name}}.broad.duplicate_line.txt

			fi
		""")

