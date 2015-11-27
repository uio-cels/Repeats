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

export REPO=/work/users/willibr/testing/repository
GENOME=$1

module load repeatmodeler/1.0.8
module load python2/2.7.9

#-----------------------------------RepeatModeler----------------------------------------#

#Making sure RepeatModeler chooses to sample from contigs larger than N50.

assemblathon_stats.pl -csv $GENOME | grep "N50 scaffold length" | column -t | \
awk '{print $4}' > N50

python2 $REPO/N50.py -i $GENOME --N50 N50 > $GENOME.N50

#Create index for RepeatModeler

BuildDatabase -name $GENOME.N50 -engine ncbi $GENOME.N50

#Run RepeatModeler

RepeatModeler -database $GENOME.N50 -engine ncbi 

