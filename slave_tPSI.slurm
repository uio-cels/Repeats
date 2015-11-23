#!/bin/sh
#SBATCH --job-name=tPSI
#SBATCH --account=nn9244k
#SBATCH --time=0:10:0
#SBATCH --mem-per-cpu=3000M 
#SBATCH --cpus-per-task=1

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi
# Mail at beginning/end/on suspension
#SBATCH --mail-user=williambrynildsen@gmail.com
#SBATCH --mail-type=ALL

GENOME=$1

#------------------------------------TransposonPSI---------------------------------------#

perl $tPSI $GENOME nuc >tPSI.log
