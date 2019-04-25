#! /bin/bash
set -o pipefail
set -e

echo "Pipeline execution initiated at: "$(date)
mkdir -p ./Script
mkdir -p ./logs/pre_process
mkdir -p ./logs/alignment
mkdir -p ./logs/post_alignment
mkdir -p ./logs/peak_calling
mkdir -p ./logs/peak_analysis
mkdir -p ./logs/peak_annotate
mkdir -p ./logs/build_trackhub

PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))

module load snakemake samtools || exit 1
module load java

snakemake --snakefile pre_process.py --configfile Encode.json --unlock
#
snakemake --snakefile pre_process.py --configfile Encode.json --cluster-config pre_process.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile alignment.py --configfile Encode.json --cluster-config alignment.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile post_alignment.py --configfile Encode.json --cluster-config post_alignment.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile peak_calling.py --configfile Encode.json --cluster-config peak_calling.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile peak_analysis.py --configfile Encode.json --cluster-config peak_analysis.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile peak_annotate.py --configfile Encode.json --cluster-config peak_annotate.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile build_trackhub.py --configfile Encode.json --cluster-config build_trackhub.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
rm -rf ./*.pos
rm -rf ./*.tmp
echo "Pipeline execution successfully finished at: "$(date)