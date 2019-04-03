# ################################### INFO ##################################### #
# Author: Amir Shams
# Date: May-15-2017
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow to process ChIp-seq data.
# ################################### IMPORT ##################################### #
import os
from os.path import join
import sys
import multiprocessing

configfile:	"config.yaml"

# LOAD INPUT/OUTPUT
INPUT_DIR = config['General_Parameters']['Input_Directory']
OUTPUT_DIR = config['General_Parameters']['Output_Directory']
EXEC_DIR = config['General_Parameters']['Exec_Directory']

# Load DATA
Report = expand(OUTPUT_DIR + "{report}_snakemake.txt", report=['Report'])

Control_List = config['FASTQ']['CONTROL']
Case_List = config['FASTQ']['CASE']
Fastq_List = Control_List + Case_List
Fastq = expand(INPUT_DIR + "{sample}.fastq", sample=Fastq_List)

Read_Alignment_DIR = OUTPUT_DIR + 'READ_ALIGNMENT/'
Read_Alignment = expand(Read_Alignment_DIR + "{sample}.bam", sample=Fastq_List)

SPP_DIR = OUTPUT_DIR + 'SPP_PEAKS/'
SPP_Peaks = expand(SPP_DIR + "{case}_vs_{control}.narrowPeak", zip, case=Case_List, control=Control_List)
SPP_IDR = expand(SPP_DIR + "{IDR}-overlapped-peaks.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])
SPP_HOMER = expand(SPP_DIR + "{IDR}_Homer_Annotate.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])

MACS_DIR = OUTPUT_DIR + 'MACS_PEAKS/'
MACS_Peaks = expand(MACS_DIR + "{case}_vs_{control}_peaks.narrowPeak", zip, case=Case_List, control=Control_List)
MACS_IDR = expand(MACS_DIR + "{IDR}-overlapped-peaks.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])
MACS_HOMER = expand(MACS_DIR + "{IDR}_Homer_Annotate.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])

GEM_DIR = OUTPUT_DIR + 'GEM_PEAKS/'
GEM_Peaks = expand(GEM_DIR + "{case}_vs_{control}.GEM_events.narrowPeak", zip, case=Case_List, control=Control_List)
GEM_IDR = expand(GEM_DIR + "{IDR}-overlapped-peaks.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])
GEM_HOMER = expand(GEM_DIR + "{IDR}_Homer_Annotate.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])

PEAKSEQ_DIR = OUTPUT_DIR + 'PEAKSEQ_PEAKS/'
PEAKSEQ_Peaks = expand(PEAKSEQ_DIR + "{case}_vs_{control}_narrowPeak.txt", zip, case=Case_List, control=Control_List)
PEAKSEQ_IDR = expand(PEAKSEQ_DIR + "{IDR}-overlapped-peaks.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])
PEAKSEQ_HOMER = expand(PEAKSEQ_DIR + "{IDR}_Homer_Annotate.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])

HOTSPOT_DIR = OUTPUT_DIR + 'HOTSPOT_PEAKS/'
HOTSPOT_Peaks = expand(HOTSPOT_DIR + "{case}_vs_{control}.hot.pks.bed", zip, case=Case_List, control=Control_List)
HOTSPOT_IDR = expand(HOTSPOT_DIR + "{IDR}-overlapped-peaks.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])
HOTSPOT_HOMER = expand(HOTSPOT_DIR + "{IDR}_Homer_Annotate.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])

JAMM_DIR = OUTPUT_DIR + 'JAMM_PEAKS/'
JAMM_Peaks = expand(JAMM_DIR + "{case}_vs_{control}.GEM_events.narrowPeak", zip, case=Case_List, control=Control_List)
JAMM_IDR = expand(JAMM_DIR + "{IDR}-overlapped-peaks.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])
JAMM_HOMER = expand(JAMM_DIR + "{IDR}_Homer_Annotate.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])

PEPR_DIR = OUTPUT_DIR + 'PEPR_PEAKS/'
PEPR_Peaks = expand(PEPR_DIR + "{case}_vs_{control}_PePr_peaks.bed", zip, case=Case_List, control=Control_List)
PEPR_IDR = expand(PEPR_DIR + "{IDR}-overlapped-peaks.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])
PEPR_HOMER = expand(PEPR_DIR + "{IDR}_Homer_Annotate.txt", IDR=['__VS__'.join(expand("{case}_vs_{control}", zip, case=Case_List, control=Control_List))])

Threads_Max = multiprocessing.cpu_count()
# ################################### EXECUTE THE WORKFLOW ####################### #


rule all:
	input: Report + Read_Alignment + SPP_Peaks + SPP_IDR + SPP_HOMER + MACS_Peaks + MACS_IDR + MACS_HOMER + GEM_Peaks + GEM_IDR + GEM_HOMER + PEAKSEQ_Peaks + PEAKSEQ_IDR + PEAKSEQ_HOMER + HOTSPOT_Peaks + HOTSPOT_IDR + HOTSPOT_HOMER + PEPR_Peaks + PEPR_IDR + PEPR_HOMER + JAMM_Peaks + JAMM_IDR + JAMM_HOMER

rule Platform_Setup:
	output:
		touch(OUTPUT_DIR + "{report}_snakemake.txt")
	priority: 10
	threads: Threads_Max
	message: "Install Required and missing Tools and packages"
	run:
		shell("""
		mkdir -p {Read_Alignment_DIR}
		mkdir -p {EXEC_DIR}
		cd {EXEC_DIR}
		echo "INSTALLING REQUIREMENTS..."
		wget --quiet http://groups.csail.mit.edu/cgs/gem/download/gem.v3.0.tar.gz -O gem.v3.0.tar.gz
		tar zxf gem.v3.0.tar.gz
		rm -rf gem.v3.0.tar.gz
		wget --quiet http://www.uwencode.org/proj/hotspot/hotspot-distr-v4.tgz -O hotspot-distr-v4.tgz
		tar zxf hotspot-distr-v4.tgz
		rm -rf hotspot-distr-v4.tgz
		mv hotspot-distr-v4 hotspot
		wget --quiet http://archive.gersteinlab.org/proj/PeakSeq/Scoring_ChIPSeq/Code/C/PeakSeq_1.31.zip -O PeakSeq_1.31.zip
		unzip -q PeakSeq_1.31.zip
		rm -rf PeakSeq_1.31.zip
		cd PeakSeq
		make --quiet
		echo "Done!"
		""")

rule Read_Alignment:
	input:
		INPUT_DIR + "{sample}.fastq"
	output:
		Read_Alignment_DIR + "{sample}.sai",
		Read_Alignment_DIR + "{sample}.sam",
		Read_Alignment_DIR + "{sample}.bam",
		Read_Alignment_DIR + "{sample}.bed",
	priority: 6
	threads: Threads_Max
	message: "Align {wildcards.sample}.clean.fastq against " + config['READ_ALIGNMENT']['REFERENCE']
	log:
		bwa_aln = Read_Alignment_DIR + "{sample}_bwa_aln_log.txt",
		bwa_samse = Read_Alignment_DIR + "{sample}_bwa_samse_log.txt",
		samtools = Read_Alignment_DIR + "{sample}_samtools_log.txt"
	run:
		if (config['READ_ALIGNMENT']['ALIGNMENT_ENGINE'] == 'BOWTIE2'):
			shell("""
			module load bowtie/2-2.2.9
			bowtie2 -x {config[READ_ALIGNMENT][BOWTIE2][BOWTIE2_INDEX} --threads={threads} {config[READ_ALIGNMENT][BOWTIE2][PRESETS]} {config[READ_ALIGNMENT][BOWTIE2][OUTPUT]} -q {input} -S {output}
			""")
		elif (config['READ_ALIGNMENT']['ALIGNMENT_ENGINE'] == 'BWA'):
			shell("""
			# ========================================
			# Map reads to create raw SAM file
			# ========================================
			module load bwa
			bwa aln -t {threads} {config[READ_ALIGNMENT][BWA][BWA_INDEX]} {input} > {output[0]} 2> {log.bwa_aln}
			bwa samse {config[READ_ALIGNMENT][BWA][BWA_INDEX]} {output[0]} {input} > {output[1]} 2> {log.bwa_samse}
			module load samtools
			samtools view --threads {threads} -Su {output[1]} -o {output[2]} > {log.samtools}
			module load bedtools
			bedtools bamtobed -i {output[2]} > {output[3]}
			""")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++SPP PEAK CALLING, IDR and ANNOTATE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule SPP_PREPROCESS:
	input:
		Case = Read_Alignment_DIR + "{case}.bam",
		Control = Read_Alignment_DIR + "{control}.bam"
	output:
		SPP_DIR + "{case}_vs_{control}_SPP_script.R"
	priority: 0
	threads: Threads_Max
	log:
		SPP_DIR + "{case}_vs_{control}_SPP_PREPROCESS_log.txt"
	message: "SPP PEAK CALLING"
	run:
		shell("""
		mkdir -p {SPP_DIR}
		echo 'library(spp);
library(snow);
Processors <- makeCluster({threads},type="SOCK");
chip.data <- read.bam.tags("{input.Case}");
input.data <- read.bam.tags("{input.Control}");
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5,cluster=Processors);
print(paste("binding peak separation distance =",binding.characteristics$peak$x))
chip.data <- select.informative.tags(chip.data,binding.characteristics);
input.data <- select.informative.tags(input.data,binding.characteristics);
chip.data <- remove.local.tag.anomalies(chip.data);
input.data <- remove.local.tag.anomalies(input.data);
fdr <- {config[SPP][FDR_RATE]};
detection.window.halfsize <- binding.characteristics$whs;
bp <- find.binding.positions(signal.data=chip.data,control.data=input.data,fdr=fdr,whs=detection.window.halfsize,cluster=Processors);
print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"));
output.binding.results(bp,"{SPP_DIR}{wildcards.case}_vs_{wildcards.control}.binding.positions.txt");
bp <- find.binding.positions(signal.data=chip.data,control.data=input.data,fdr=fdr,method=tag.lwcc,whs=detection.window.halfsize,cluster=Processors);
bp <- add.broad.peak.regions(chip.data,input.data,bp,window.size=1000,z.thr=3);
write.narrowpeak.binding(bp,"{SPP_DIR}{wildcards.case}_vs_{wildcards.control}.narrowPeak")' > {output}
""")


rule SPP_PEAK_CALLING:
	input:
		SPP_DIR + "{case}_vs_{control}_SPP_script.R"
	output:
		SPP_DIR + "{case}_vs_{control}.narrowPeak",
		SPP_DIR + "{case}_vs_{control}.binding.positions.txt"
	priority: 0
	threads: Threads_Max
	message: "SPP PEAK CALLING"
	run:
		shell("""
			module load R
			module load samtools
			Rscript {input}
			""")
rule SPP_IDR:
	input:
		expand(SPP_DIR + "{case}_vs_{control}.narrowPeak", zip, case=Case_List, control=Control_List)
	output:
		SPP_DIR + "{IDR}-overlapped-peaks.txt"
	priority: 0
	threads: Threads_Max
	message: "SPP IDR"
	log:
		SPP_DIR + "{IDR}-IDR_Report.txt"
	run:
		shell("""
			module load idr
			idr --samples {input} --input-file-type {config[SPP][IDR][PEAK_TYPE]} --output-file {output[0]} --rank {config[SPP][IDR][RANK]} --soft-idr-threshold {config[SPP][IDR][SOFT_IDR_THRESHOLD]} \
--peak-merge-method {config[SPP][IDR][PEAK_MERGE_METHOD]} --initial-mu {config[SPP][IDR][INITIAL_MU]} --initial-sigma {config[SPP][IDR][INITIAL_SIGMA]} --initial-rho {config[SPP][IDR][INITIAL_RHO]} \
--initial-mix-param  {config[SPP][IDR][INITIAL_MIX_PARAM]} --random-seed {config[SPP][IDR][RANDOM_SEED]} --max-iter {config[SPP][IDR][MAX_ITER]} --convergence-eps {config[SPP][IDR][CONVERGENCE_EPS]} \
{config[SPP][IDR][BEST_MULTI_SUMMIT]} {config[SPP][IDR][VERBOSE]} --plot --log-output-file {log}
			""")

rule SPP_HOMER:
	input:
		SPP_DIR + "{IDR}-overlapped-peaks.txt"
	output:
		SPP_DIR + "{IDR}_Homer_Annotate.txt"
	log:
		SPP_DIR + "{IDR}_Homer_Annotate_stats.txt"
	priority: 0
	threads: Threads_Max
	message: "SPP PEAK ANNOTATION"
	run:
		shell("""
			module load homer
			annotatePeaks.pl {input} hg19 > {output} 2>{log}
			""")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++MACS PEAK CALLING, IDR and ANNOTATE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule MACS_PEAKS_CALLING:
	input:
		Case = Read_Alignment_DIR + "{case}.bam",
		Control = Read_Alignment_DIR + "{control}.bam"
	output:
		MACS_DIR + "{case}_vs_{control}_model.r",
		MACS_DIR + "{case}_vs_{control}_peaks.narrowPeak",
		MACS_DIR + "{case}_vs_{control}_peaks.xls",
		MACS_DIR + "{case}_vs_{control}_summits.bed",
		MACS_DIR + "{case}_vs_{control}_model.pdf"
	log:
		MACS_DIR + "{case}_vs_{control}_MACS_CALL_PEAK_REPORT.txt"
	priority: 3
	threads: Threads_Max
	message: "MACS PEAK CALLINGS"
	run:
		shell("""
			mkdir -p {MACS_DIR}
			module load macs/2.1.0.20150420
			macs2 callpeak --treatment {input.Case} --control {input.Control} --format {config[MACS][FORMAT]}  --gsize {config[MACS][G_SIZE]}  --keep-dup {config[MACS][KEEP_DUPLICATES]} --buffer-size {config[MACS][BUFFER_SIZE]} \
--outdir {MACS_DIR} --name '{wildcards.case}_vs_{wildcards.control}' --bw {config[MACS][BANDWIDTH]} --verbose {config[MACS][VERBOSE]} --shift {config[MACS][SHIFT]} --extsize {config[MACS][EXTSIZE]} --qvalue {config[MACS][Q_VALUE_CUTOFF]} \
--fe-cutoff {config[MACS][FOLD_ENRICHMENT_CUTOFF]}  --slocal {config[MACS][SMALL_LOCAL]} --llocal {config[MACS][LARGE_LOCAL]} --mfold {config[MACS][M_FOLD]} {config[MACS][TO_LARGE]} 2> {log}
			module load R
			cd {MACS_DIR}
			Rscript {output[0]}
			""")

rule MACS_IDR:
	input:
		expand(MACS_DIR + "{case}_vs_{control}_peaks.narrowPeak", zip, case=Case_List, control=Control_List)
	output:
		MACS_DIR + "{IDR}-overlapped-peaks.txt"
	priority: 0
	threads: Threads_Max
	message: "MACS IDR"
	log:
		MACS_DIR + "{IDR}-IDR_Report.txt"
	run:
		shell("""
			module load idr
			idr --samples {input} --input-file-type {config[MACS][IDR][PEAK_TYPE]} --output-file {output[0]} --rank {config[MACS][IDR][RANK]} --soft-idr-threshold {config[MACS][IDR][SOFT_IDR_THRESHOLD]} \
--peak-merge-method {config[MACS][IDR][PEAK_MERGE_METHOD]} --initial-mu {config[MACS][IDR][INITIAL_MU]} --initial-sigma {config[MACS][IDR][INITIAL_SIGMA]} --initial-rho {config[MACS][IDR][INITIAL_RHO]} \
--initial-mix-param  {config[MACS][IDR][INITIAL_MIX_PARAM]} --random-seed {config[MACS][IDR][RANDOM_SEED]} --max-iter {config[MACS][IDR][MAX_ITER]} --convergence-eps {config[MACS][IDR][CONVERGENCE_EPS]} \
{config[MACS][IDR][BEST_MULTI_SUMMIT]} {config[MACS][IDR][VERBOSE]} --plot --log-output-file {log}
			""")

rule MACS_HOMER:
	input:
		MACS_DIR + "{IDR}-overlapped-peaks.txt"
	output:
		MACS_DIR + "{IDR}_Homer_Annotate.txt"
	log:
		MACS_DIR + "{IDR}_Homer_Annotate_stats.txt"
	priority: 0
	threads: Threads_Max
	message: "MACS PEAK ANNOTATION"
	run:
		shell("""
			module load homer
			annotatePeaks.pl {input} hg19 > {output} 2>{log}
			""")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++GEM PEAK CALLING, IDR and ANNOTATE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule GEM_PEAKS_CALLING:
	input:
		Case = Read_Alignment_DIR + "{case}.bed",
		Control = Read_Alignment_DIR + "{control}.bed"
	output:
		GEM_DIR + "{case}_vs_{control}.GEM_events.narrowPeak"
	threads: Threads_Max
	log:
		GEM_DIR + "{case}_vs_{control}_GEM_CALL_PEAK_REPORT.txt"
	message: "GEM PEAK CALLING"
	run:
		shell("""
			# ================================
			# GEM
			# ================================
			mkdir -p {GEM_DIR}
			cd {EXEC_DIR}gem
			nohup java -XX:ParallelGCThreads={threads} -jar gem.jar --d Read_Distribution_default.txt --g hg19.chrom.sizes  --s {config[GEM][S]} --expt {input.Case}  --ctrl {input.Control} --f {config[GEM][F]} --out {GEM_DIR}{wildcards.case}_vs_{wildcards.control} --genome {config[GEM][GENOME]} --k_min {config[GEM][KMIN]} --k_max {config[GEM][KMAX]} --outNP --q {config[GEM][Q]}
			sort -k7nr,7nr {GEM_DIR}{wildcards.case}_vs_{wildcards.control}/{wildcards.case}_vs_{wildcards.control}.GEM_events.narrowPeak  > {output} &
			mv nohup {log}
			""")

rule GEM_IDR:
	input:
		expand(GEM_DIR + "{case}_vs_{control}.GEM_events.narrowPeak", zip, case=Case_List, control=Control_List)
	output:
		GEM_DIR + "{IDR}-overlapped-peaks.txt"
	priority: 0
	threads: Threads_Max
	message: "GEM IDR"
	log:
		GEM_DIR + "{IDR}-IDR_Report.txt"
	run:
		shell("""
			module load idr
			idr --samples {input} --input-file-type {config[GEM][IDR][PEAK_TYPE]} --output-file {output[0]} --rank {config[GEM][IDR][RANK]} --soft-idr-threshold {config[GEM][IDR][SOFT_IDR_THRESHOLD]} \
--peak-merge-method {config[GEM][IDR][PEAK_MERGE_METHOD]} --initial-mu {config[GEM][IDR][INITIAL_MU]} --initial-sigma {config[GEM][IDR][INITIAL_SIGMA]} --initial-rho {config[GEM][IDR][INITIAL_RHO]} \
--initial-mix-param  {config[GEM][IDR][INITIAL_MIX_PARAM]} --random-seed {config[GEM][IDR][RANDOM_SEED]} --max-iter {config[GEM][IDR][MAX_ITER]} --convergence-eps {config[GEM][IDR][CONVERGENCE_EPS]} \
{config[GEM][IDR][BEST_MULTI_SUMMIT]} {config[GEM][IDR][VERBOSE]} --plot --log-output-file {log}
			""")

rule GEM_HOMER:
	input:
		GEM_DIR + "{IDR}-overlapped-peaks.txt"
	output:
		GEM_DIR + "{IDR}_Homer_Annotate.txt"
	log:
		GEM_DIR + "{IDR}_Homer_Annotate_stats.txt"
	priority: 0
	threads: Threads_Max
	message: "GEM PEAK ANNOTATION"
	run:
		shell("""
			module load homer
			annotatePeaks.pl {input} hg19 > {output} 2>{log}
			""")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++PEAKSEQ PEAK CALLING, IDR and ANNOTATE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule PEAKSEQ_PREPROCESS:
	input:
		Case = Read_Alignment_DIR + "{case}.bed",
		Control = Read_Alignment_DIR + "{control}.bed"
	output:
		PEAKSEQ_DIR + "{case}_vs_{control}_config.dat"
	threads: Threads_Max
	message: "PEAKSEQ PREPROCESS"
	run:
		shell("""
			# ================================
			# PEAKSEQ PREPROCESS
			# ================================
			mkdir -p {PEAKSEQ_DIR}
			mkdir -p {PEAKSEQ_DIR}Peak_Seq_{wildcards.case}
			mkdir -p {PEAKSEQ_DIR}Peak_Seq_{wildcards.control}
			{EXEC_DIR}PeakSeq/bin/PeakSeq -preprocess tagAlign {input.Case} {PEAKSEQ_DIR}Peak_Seq_{wildcards.case}
			{EXEC_DIR}PeakSeq/bin/PeakSeq -preprocess tagAlign {input.Control} {PEAKSEQ_DIR}Peak_Seq_{wildcards.control}
			# ==========================================================================
			echo "# Experiment id is used as a prefix to the output file name.
Experiment_id {wildcards.case}_vs_{wildcards.control}
# Enrichment fragment length For tag extension, this is the value of average fragment length.
Enrichment_mapped_fragment_length {config[PEAKSEQ][FRAGMENT_LENGTH]}
# Target FDR in the simulations.
target_FDR {config[PEAKSEQ][TARGET_FDR]}
# Number of simulations performed while estimating the putative peaks.
N_Simulations {config[PEAKSEQ][N_SIMULATION]}
# Minimum distance between consecutive peaks
Minimum_interpeak_distance {config[PEAKSEQ][MIN_INTERPEAK_DISTANCE]}
#The directory that contains the preprocessed ChIP-Seq reads, can specify multiple directories to pool reads from multiple source (e.g. replicates)
ChIP_Seq_reads_data_dirs {PEAKSEQ_DIR}Peak_Seq_{wildcards.case}
#The directory that contains the preprocessed Input (control) experiment reads. (Multiple directories allowed)
Input_reads_data_dirs {PEAKSEQ_DIR}Peak_Seq_{wildcards.control}
# Mappability file that includes the uniquely mappable number of nucleotides per window for each chromosome.
Mappability_map_file {config[PEAKSEQ][MAP]}
#Seed for pseudo-random number generator. This is necessary for simulated background option (specified below).
#Simulation_seed 1234567
#Q-value threshold applied on the final set of peaks.
max_Qvalue {config[PEAKSEQ][MAX_QVALUE]}
# There are currently two models for simulating the background for threshold selection
# Simulated background is the simulation based method that is explained in the PeakSeq paper.
# Poisson background uses a simple Poisson background with mean estimated from the read statistics. This option is still experimental but it is much faster than the simulated background option.
# Background_model Poisson
Background_model {config[PEAKSEQ][BACKGROUND_MODEL]}" > {output}
""")

rule PEAKSEQ_PEAK_CALLING:
	input:
		PEAKSEQ_DIR + "{case}_vs_{control}_config.dat"
	output:
		PEAKSEQ_DIR + "{case}_vs_{control}_narrowPeak.txt"
	threads: Threads_Max
	message: "PEAKSEQ PEAK CALLING"
	log:
		PEAKSEQ_DIR + "{case}_vs_{control}_PEAKSEQ_log.txt"
	run:
		shell("""
			cd {PEAKSEQ_DIR}
			{EXEC_DIR}PeakSeq/bin/PeakSeq -peak_select {input} > {log} 2>&1
			""")

rule PEAKSEQ_IDR:
	input:
		expand(PEAKSEQ_DIR + "{case}_vs_{control}_narrowPeak.txt", zip, case=Case_List, control=Control_List)
	output:
		PEAKSEQ_DIR + "{IDR}-overlapped-peaks.txt"
	threads: Threads_Max
	message: "PEAKSEQ IDR"
	log:
		PEAKSEQ_DIR + "{IDR}-IDR_Report.txt"
	run:
		shell("""
			module load idr
			idr --samples {input} --input-file-type {config[PEAKSEQ][IDR][PEAK_TYPE]} --output-file {output[0]} --rank {config[PEAKSEQ][IDR][RANK]} --soft-idr-threshold {config[PEAKSEQ][IDR][SOFT_IDR_THRESHOLD]} \
--peak-merge-method {config[PEAKSEQ][IDR][PEAK_MERGE_METHOD]} --initial-mu {config[PEAKSEQ][IDR][INITIAL_MU]} --initial-sigma {config[PEAKSEQ][IDR][INITIAL_SIGMA]} --initial-rho {config[PEAKSEQ][IDR][INITIAL_RHO]} \
--initial-mix-param  {config[PEAKSEQ][IDR][INITIAL_MIX_PARAM]} --random-seed {config[PEAKSEQ][IDR][RANDOM_SEED]} --max-iter {config[PEAKSEQ][IDR][MAX_ITER]} --convergence-eps {config[PEAKSEQ][IDR][CONVERGENCE_EPS]} \
{config[PEAKSEQ][IDR][BEST_MULTI_SUMMIT]} {config[PEAKSEQ][IDR][VERBOSE]} --plot --log-output-file {log}
			""")

rule PEAKSEQ_HOMER:
	input:
		PEAKSEQ_DIR + "{IDR}-overlapped-peaks.txt"
	output:
		PEAKSEQ_DIR + "{IDR}_Homer_Annotate.txt"
	log:
		PEAKSEQ_DIR + "{IDR}_Homer_Annotate_stats.txt"
	threads: Threads_Max
	message: "PEAKSEQ PEAK ANNOTATION"
	run:
		shell("""
			module load homer
			annotatePeaks.pl {input} hg19 > {output} 2> {log}
			""")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++HOTSPOT PEAK CALLING, IDR and ANNOTATE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule HOTSPOT_PRE_PROCESSING:
	input:
		Case = Read_Alignment_DIR + "{case}.bam",
		Control = Read_Alignment_DIR + "{control}.bam"
	output:
		HOTSPOT_DIR + "{case}_vs_{control}/hotspot_config.txt",
		HOTSPOT_DIR + "{case}_vs_{control}/hotspot_peakcalling.sh"
	threads: Threads_Max
	message: "HOTSPOT PREPROCESS"
	run:
		shell("""
			# ================================
			# HOTSPOT
			# ================================
			mkdir -p {HOTSPOT_DIR}
			cd {HOTSPOT_DIR}
			mkdir -p {wildcards.case}_vs_{wildcards.control}
			# ==========================================================================
			echo "[script-tokenizer]
# HOTSPOT PARAMETERS FILE:
_TAGS_ = {input.Case}
_USE_INPUT_ = F
_INPUT_TAGS_ =
_GENOME_ = {config[HOTSPOT][GENOME]}
_K_ = {config[HOTSPOT][TAG_LENGTH]}
_CHROM_FILE_ = {EXEC_DIR}hotspot/data/hg19.chromInfo.bed
_MAPPABLE_FILE_ = {EXEC_DIR}hotspot/data/hg19.K36.mappable_only.bed.starch
_DUPOK_ = {config[HOTSPOT][DUP_OK]}
_FDRS_ = {config[HOTSPOT][FDRS]}
_DENS_:
_OUTDIR_ = {HOTSPOT_DIR}/{wildcards.case}_vs_{wildcards.control}
_RANDIR_ = {HOTSPOT_DIR}/{wildcards.case}_vs_{wildcards.control}
_OMIT_REGIONS_: {EXEC_DIR}hotspot/data/Satellite.hg19.bed
_CHECK_ = F
_CHKCHR_ = chrX
_HOTSPOT_ = {EXEC_DIR}hotspot/hotspot-deploy/bin/hotspot
_CLEAN_ = T
_PKFIND_BIN_ = {EXEC_DIR}hotspot/hotspot-deploy/bin/wavePeaks
_PKFIND_SMTH_LVL_ = {config[HOTSPOT][PKFIND_SMTH_LVL]}
_SEED_=101
## Hotspot program parameters
_THRESH_ = {config[HOTSPOT][THRESH]}
_WIN_MIN_ = {config[HOTSPOT][WIN_MIN]}
_WIN_MAX_ = {config[HOTSPOT][WIN_MAX]}
_WIN_INCR_ = {config[HOTSPOT][WIN_INCR]}
_BACKGRD_WIN_ = {config[HOTSPOT][BACKGRD_WIN]}
_MERGE_DIST_ = {config[HOTSPOT][MERGE_DIST]}
_MINSIZE_ = {config[HOTSPOT][MINSIZE]}" > {output[0]}
echo "###################################################"
echo '#! /bin/bash
scriptTokBin={EXEC_DIR}hotspot/ScriptTokenizer/src/script-tokenizer.py
pipeDir={EXEC_DIR}hotspot/pipeline-scripts
tokenFile={output[0]}
## Do SPOT only (set _FDRS_ to "N" in runall.tokens.txt)
# scripts="$pipeDir/run_make_lib
#     $pipeDir/run_10kb_counts
#     $pipeDir/run_pass1_hotspot
#     $pipeDir/run_pass1_merge_and_thresh_hotspots
#     $pipeDir/run_pass2_hotspot
#     $pipeDir/run_rescore_hotspot_passes
#     $pipeDir/run_spot"
## Do everything, including badspots and final cleanup
scripts="$pipeDir/run_badspot
	$pipeDir/run_make_lib
	$pipeDir/run_wavelet_peak_finding
	$pipeDir/run_10kb_counts
	$pipeDir/run_generate_random_lib
	$pipeDir/run_pass1_hotspot
	$pipeDir/run_pass1_merge_and_thresh_hotspots
	$pipeDir/run_pass2_hotspot
	$pipeDir/run_rescore_hotspot_passes
	$pipeDir/run_spot
	$pipeDir/run_thresh_hot.R
	$pipeDir/run_both-passes_merge_and_thresh_hotspots
	$pipeDir/run_add_peaks_per_hotspot
	$pipeDir/run_final"
$scriptTokBin \
	--clobber \
	--output-dir=`pwd` \
	$tokenFile \
	$scripts
module load hotspot bedops bedtools R
module load R
for script in $scripts
do
	./$(basename $script).tok
done' > {output[1]}
chmod 755 {output[1]}
""")

rule HOTSPOT_PEAK_CALLING:
	input:
		HOTSPOT_DIR + "{case}_vs_{control}/hotspot_peakcalling.sh"
	output:
		HOTSPOT_DIR + "{case}_vs_{control}/{case}-final/{case}.hot.bed"
	threads: Threads_Max
	message: "HOTSPOT PEAK CALLING"
	log:
		HOTSPOT_DIR + "{case}_vs_{control}_HOTSPOT_log.txt"
	run:
		shell("""
			cd {HOTSPOT_DIR}{wildcards.case}_vs_{wildcards.control}
			bash ./hotspot_peakcalling.sh > {log} 2>&1
			""")

rule HOTSPOT_POST_PROCESS:
	input:
		HOTSPOT_DIR + "{case}_vs_{control}/{case}-final/{case}.hot.bed"
	output:
		HOTSPOT_DIR + "{case}_vs_{control}.hot.pks.bed"
	threads: Threads_Max
	message: "HOTSPOT IDR"
	run:
		shell("""
			cd {HOTSPOT_DIR}{wildcards.case}_vs_{wildcards.control}
			cd {wildcards.case}-final
			sed \"s/$/\t\./\" {wildcards.case}.hot.bed > {wildcards.case}.hot.pks.bed
			paste {wildcards.case}.hot.pks.bed {wildcards.case}.hot.pval.txt > temp.txt
			sed \"s/$/\t\-1\t-1/\" temp.txt > {wildcards.case}.hot.pks.bed
			cp {wildcards.case}.hot.pks.bed {HOTSPOT_DIR}{wildcards.case}_vs_{wildcards.control}.hot.pks.bed
			""")

rule HOTSPOT_IDR:
	input:
		expand(HOTSPOT_DIR + "{case}_vs_{control}.hot.pks.bed", zip, case=Case_List, control=Control_List)
	output:
		HOTSPOT_DIR + "{IDR}-overlapped-peaks.txt"
	threads: Threads_Max
	message: "HOTSPOT IDR"
	log:
		HOTSPOT_DIR + "{IDR}-IDR_Report.txt"
	run:
		shell("""
			module load idr
			idr --samples {input} --input-file-type {config[HOTSPOT][IDR][PEAK_TYPE]} --output-file {output[0]} --rank {config[HOTSPOT][IDR][RANK]} --soft-idr-threshold {config[HOTSPOT][IDR][SOFT_IDR_THRESHOLD]} \
--peak-merge-method {config[HOTSPOT][IDR][PEAK_MERGE_METHOD]} --initial-mu {config[HOTSPOT][IDR][INITIAL_MU]} --initial-sigma {config[HOTSPOT][IDR][INITIAL_SIGMA]} --initial-rho {config[HOTSPOT][IDR][INITIAL_RHO]} \
--initial-mix-param  {config[HOTSPOT][IDR][INITIAL_MIX_PARAM]} --random-seed {config[HOTSPOT][IDR][RANDOM_SEED]} --max-iter {config[HOTSPOT][IDR][MAX_ITER]} --convergence-eps {config[HOTSPOT][IDR][CONVERGENCE_EPS]} \
{config[HOTSPOT][IDR][BEST_MULTI_SUMMIT]} {config[HOTSPOT][IDR][VERBOSE]} --plot --log-output-file {log}
			""")

rule HOTSPOT_HOMER:
	input:
		HOTSPOT_DIR + "{IDR}-overlapped-peaks.txt"
	output:
		HOTSPOT_DIR + "{IDR}_Homer_Annotate.txt"
	log:
		HOTSPOT_DIR + "{IDR}_Homer_Annotate_stats.txt"
	threads: Threads_Max
	message: "HOTSPOT PEAK ANNOTATION"
	run:
		shell("""
			module load homer
			annotatePeaks.pl {input} hg19 > {output} 2>{log}
			""")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++JAMM PEAK CALLING, IDR and ANNOTATE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule JAMM_PEAKS_CALLING:
	input:
		Case = Read_Alignment_DIR + "{case}.bed",
		Control = Read_Alignment_DIR + "{control}.bed"
	output:
		JAMM_DIR + "{case}_vs_{control}.GEM_events.narrowPeak"
	threads: Threads_Max
	log:
		JAMM_DIR + "{case}_vs_{control}_JAMM_CALL_PEAK_REPORT.txt"
	message: "JAMM PEAK CALLING"
	run:
		shell("""
			# ================================
			# JAMM
			# ================================
			mkdir -p {JAMM_DIR}
			cd {JAMM_DIR}
			mkdir -p {wildcards.case}_vs_{wildcards.control}
			cd {wildcards.case}_vs_{wildcards.control}
			mkdir -p {wildcards.case}
			mkdir -p {wildcards.control}
			cp {input.Case} {wildcards.case}
			cp {input.Control} {wildcards.control}
			wget --quiet https://raw.githubusercontent.com/amirshams84/Chips/master/jamm_hg19.txt
			module load JAMM
			JAMM.sh -s ./{wildcards.case} -c ./{wildcards.control} -g ./jamm_hg19.txt -o ./{wildcards.case}_vs_{wildcards.control} -r peak -m normal -p {Threads_Max} > {log} 2>&1
			""")

rule JAMM_IDR:
	input:
		expand(JAMM_DIR + "{case}_vs_{control}.GEM_events.narrowPeak", zip, case=Case_List, control=Control_List)
	output:
		JAMM_DIR + "{IDR}-overlapped-peaks.txt"
	threads: Threads_Max
	message: "JAMM IDR"
	log:
		JAMM_DIR + "{IDR}-IDR_Report.txt"
	run:
		shell("""
			module load idr
			idr --samples {input} --input-file-type {config[JAMM][IDR][PEAK_TYPE]} --output-file {output[0]} --rank {config[JAMM][IDR][RANK]} --soft-idr-threshold {config[JAMM][IDR][SOFT_IDR_THRESHOLD]} \
--peak-merge-method {config[JAMM][IDR][PEAK_MERGE_METHOD]} --initial-mu {config[JAMM][IDR][INITIAL_MU]} --initial-sigma {config[JAMM][IDR][INITIAL_SIGMA]} --initial-rho {config[JAMM][IDR][INITIAL_RHO]} \
--initial-mix-param  {config[JAMM][IDR][INITIAL_MIX_PARAM]} --random-seed {config[JAMM][IDR][RANDOM_SEED]} --max-iter {config[JAMM][IDR][MAX_ITER]} --convergence-eps {config[JAMM][IDR][CONVERGENCE_EPS]} \
{config[JAMM][IDR][BEST_MULTI_SUMMIT]} {config[JAMM][IDR][VERBOSE]} --plot --log-output-file {log}
			""")

rule JAMM_HOMER:
	input:
		JAMM_DIR + "{IDR}-overlapped-peaks.txt"
	output:
		JAMM_DIR + "{IDR}_Homer_Annotate.txt"
	log:
		JAMM_DIR + "{IDR}_Homer_Annotate_stats.txt"
	priority: 0
	threads: Threads_Max
	message: "JAMM PEAK ANNOTATION"
	run:
		shell("""
			module load homer
			annotatePeaks.pl {input} hg19 > {output} 2>{log}
			""")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++PEPR PEAK CALLING, IDR and ANNOTATE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rule PEPR_PRE_PROCESSING:
	input:
		Case = Read_Alignment_DIR + "{case}.bed",
		Control = Read_Alignment_DIR + "{control}.bed"
	output:
		PEPR_DIR + "{case}_vs_{control}/{case}_Replicate1.bed",
		PEPR_DIR + "{case}_vs_{control}/{case}_Replicate2.bed",
		PEPR_DIR + "{case}_vs_{control}/{control}_Replicate1.bed",
		PEPR_DIR + "{case}_vs_{control}/{control}_Replicate2.bed"
	threads: Threads_Max
	message: "PEPR PREPROCESS"
	run:
		shell("""
			# ================================
			# PEPR PREPROCESS
			# ================================
			mkdir -p {PEPR_DIR}
			cd {PEPR_DIR}
			mkdir -p {wildcards.case}_vs_{wildcards.control}
			cd {wildcards.case}_vs_{wildcards.control}
			nlines=$(cat {input.Case} | wc -l)
			nlines=$(((nlines+1)/2))
			cat {input.Case} | shuf | split -d -l $nlines - {wildcards.case}
			mv {wildcards.case}00 {output[0]}
			mv {wildcards.case}01 {output[1]}
			nlines=$(cat {input.Control} | wc -l)
			nlines=$(((nlines + 1)/2))
			cat {input.Control} | shuf | split -d -l $nlines - {wildcards.control}
			mv {wildcards.control}00 {output[2]}
			mv {wildcards.control}01 {output[3]}
			""")

rule PEPR_PEAKS_CALLING:
	input:
		PEPR_DIR + "{case}_vs_{control}/{case}_Replicate1.bed",
		PEPR_DIR + "{case}_vs_{control}/{case}_Replicate2.bed",
		PEPR_DIR + "{case}_vs_{control}/{control}_Replicate1.bed",
		PEPR_DIR + "{case}_vs_{control}/{control}_Replicate2.bed"
	output:
		PEPR_DIR + "{case}_vs_{control}_PePr_peaks.bed"
	threads: Threads_Max
	log:
		PEPR_DIR + "{case}_vs_{control}_PEPR_CALL_PEAK_REPORT.txt"
	message: "PEPR PEAK CALLING"
	run:
		shell("""
			# ================================
			# PEPR
			# ================================
			cd {PEPR_DIR}
			module load PePr
			PePr -c {input[0]},{input[1]} -i {input[2]},{input[3]} -f bed --peaktype sharp --num-processors {Threads_Max} --output-directory {wildcards.case}_vs_{wildcards.control} > {log} 2>&1
			mv {wildcards.case}_vs_{wildcards.control}/NA__PePr_peaks.bed {output}
			""")

rule PEPR_IDR:
	input:
		expand(PEPR_DIR + "{case}_vs_{control}_PePr_peaks.bed", zip, case=Case_List, control=Control_List)
	output:
		PEPR_DIR + "{IDR}-overlapped-peaks.txt"
	threads: Threads_Max
	message: "PEPR IDR"
	log:
		PEPR_DIR + "{IDR}-IDR_Report.txt"
	run:
		shell("""
			module load idr
			idr --samples {input} --input-file-type {config[PEPR][IDR][PEAK_TYPE]} --output-file {output[0]} --rank {config[PEPR][IDR][RANK]} --soft-idr-threshold {config[PEPR][IDR][SOFT_IDR_THRESHOLD]} \
--peak-merge-method {config[PEPR][IDR][PEAK_MERGE_METHOD]} --initial-mu {config[PEPR][IDR][INITIAL_MU]} --initial-sigma {config[PEPR][IDR][INITIAL_SIGMA]} --initial-rho {config[PEPR][IDR][INITIAL_RHO]} \
--initial-mix-param  {config[PEPR][IDR][INITIAL_MIX_PARAM]} --random-seed {config[PEPR][IDR][RANDOM_SEED]} --max-iter {config[PEPR][IDR][MAX_ITER]} --convergence-eps {config[PEPR][IDR][CONVERGENCE_EPS]} \
{config[PEPR][IDR][BEST_MULTI_SUMMIT]} {config[PEPR][IDR][VERBOSE]} --plot --log-output-file {log}
			""")

rule PEPR_HOMER:
	input:
		PEPR_DIR + "{IDR}-overlapped-peaks.txt"
	output:
		PEPR_DIR + "{IDR}_Homer_Annotate.txt"
	log:
		PEPR_DIR + "{IDR}_Homer_Annotate_stats.txt"
	priority: 0
	threads: Threads_Max
	message: "PEPR PEAK ANNOTATION"
	run:
		shell("""
			module load homer
			annotatePeaks.pl {input} hg19 > {output} 2>{log}
			""")
