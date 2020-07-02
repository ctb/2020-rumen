#!/bin/bash -login
#SBATCH -p med2
#SBATCH -J rumen
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org
#SBATCH -t 1-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=135gb

. "/home/ctbrown/miniconda3/etc/profile.d/conda.sh"

cd /home/ctbrown/2020-rumen

conda activate sgc-run

#spacegraphcats sgc-conf.yml search
snakemake -p --use-conda -j 4

set -o nounset
set -o errexit
set -x

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

