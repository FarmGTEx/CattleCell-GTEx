#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 16                       # Number of CPU cores
#SBATCH --mem=20480               # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH --account cattle_gtexs     #project name
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -J TWAS_model                # Name of the job
#SBATCH -t 10:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#

###Call the corresponding R version
#module load R/4.1.0-gcc-4.8.5
###Set workspace
#work_path=$1
#type=$2
#covar_path=$3
#software_path=$4
#cd ${work_path}
###Extraction string
#tissue_list=(`ls ${covar_path}`)
#Tissues=($(echo ${tissue_list[*]} | sed 's/.covariates.txt//g' | sort -u))

###Cycle submit task
#The sample size of muscle is too large. Run it at another node

base_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/"
for tissue_dir in "$base_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  for ct_dir in "$tissue_dir"/*/; do
    ct_name="$(basename "$ct_dir")"
    for chrom in {1..29}; do
      ~/miniconda3/envs/R/bin/Rscript /faststorage/project/cattle_gtexs/script/TWAS/makeModels/8.make_db.R "$tissue_name" "$ct_name"
    done
  done
done

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
