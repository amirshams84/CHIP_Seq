#! /bin/bash
set -o pipefail
set -e

echo "Pipeline execution initiated at: "$(date)

mkdir -p ./logs/pre_process
mkdir -p ./logs/alignment
mkdir -p ./logs/post_alignment
mkdir -p ./logs/peak_calling
mkdir -p ./logs/signal

PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))

module load snakemake samtools || exit 1
module load java

snakemake --snakefile PreProcess.py --configfile Yoko.json --unlock
#
snakemake --snakefile PreProcess.py --configfile Yoko.json --cluster-config PreProcess.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile Alignment.py --configfile Yoko.json --cluster-config Alignment.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile Post_Alignment.py --configfile Yoko.json --cluster-config Post_Alignment.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile Peak_Calling.py --configfile Yoko.json --cluster-config Peak_Calling.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"
#
snakemake --snakefile Signal.py --configfile Yoko.json --cluster-config Signal.yaml -j 100 --latency-wait 120 --keep-going --local-cores=$PROCESSORS --cores=$PROCESSORS --keep-going --rerun-incomplete --cluster="sbatch -c {threads} \
--mem={resources.mem_mb} --partition={cluster.partition} --job-name={cluster.jobname} --output={cluster.output} --error={cluster.error} --time={cluster.time} {cluster.extra}"


echo "Pipeline execution successfully finished at: "$(date)