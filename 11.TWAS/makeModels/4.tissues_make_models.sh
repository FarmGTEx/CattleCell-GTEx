#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
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

#base_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/TWAS/"
#for tissue_dir in "$base_dir"/*/; do
  #tissue_name="$(basename "$tissue_dir")"
  #for ct_dir in "$tissue_dir"/*/; do
    #ct_name="$(basename "$ct_dir")"
    #for chrom in {1..29}; do
      #sbatch --export=chrom="$chrom",ct_name="$ct_name",tissue_name="$tissue_name" run_TWAS_model.sh
    #done
  #done
#done

base_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
for tissue_dir in "$base_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  for chrom in {1..29}; do
    sbatch --export=chrom="$chrom",tissue_name="$tissue_name" run_TWAS_model.sh
    #~/miniconda3/envs/R/bin/Rscript /faststorage/project/cattle_gtexs/script/TWAS/makeModels/5.parallel_tiss_chrom_training.R "$chrom" "$tissue_name"
  done
done
#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID