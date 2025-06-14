#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=204800               # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH --account cattle_gtexs     #project name
#SBATCH -J plink               # Name of the job
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -t 10:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#

cd /faststorage/project/cattle_gtexs/CattleGTEx/Ancient_DNA
vcftools --gzvcf ancient_cattle_20.vcf.gz --weir-fst-pop modern_euro.txt --weir-fst-pop ancient_euro.txt --out euro_fst.txt --fst-window-size 30000 --fst-window-step 10000

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID