#!/bin/sh
#SBATCH --job-name=repmod
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

module load repeatmodeler/1.0.8

#-----------------------------------RepeatModeler----------------------------------------#

#Create index for RepeatModeler

BuildDatabase -name $GENOME -engine ncbi $GENOME

#Run RepeatModeler

RepeatModeler -database $GENOME -engine ncbi

