#!/bin/bash
set -e

module load python/3.8
module load snakemake/5.24.1

R=$2

mkdir -p ${R}/snakejobs
mkdir -p $R/Reports

if [ $1 == "npr" ]
then
    snakemake --unlock --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
    snakemake -npr --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
fi

if [ $1 == "process" ]
then
### WORKING
snakemake --unlock --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
snakemake --latency-wait 120 --configfile $R/config.yaml -s $R/Snakefile -d $R --printshellcmds --use-conda --cluster-config $R/cluster.json --keep-going --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e snakejobs/slurm-%j_{params.rname}.out -o snakejobs/slurm-%j_{params.rname}.out" -j 500 --rerun-incomplete --stats $R/Reports/snakemake.stats | tee -a $R/Reports/snakemake.log
fi
