#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 16                       # Number of CPU cores
#SBATCH --mem=20480               # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J TWAS                # Name of the job
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


base_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/"
gwas_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/gwas/"
for tissue_dir in "$base_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  for ct_dir in "$tissue_dir"/*/; do
    ct_name="$(basename "$ct_dir")"
    for trait_file in "$gwas_dir"/*; do
      trait_name="$(basename "$trait_file")"
      python3 /faststorage/project/cattle_gtexs/script/TWAS/runSoftware/SPrediXcan.py \
        --model_db_path /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/${tissue_name}/${ct_name}/dbs/eQTL_${tissue_name}_${ct_name}_ElasticNet_models_filter_0.05.db \
        --covariance /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/${tissue_name}/${ct_name}/covariances/eQTL.${tissue_name}_${ct_name}_Model_training_covariances.txt.gz \
        --gwas_file ${trait_file} \
        --snp_column SNP \
        --effect_allele_column A1 \
        --non_effect_allele_column A2 \
        --beta_column BETA \
        --pvalue_column P \
        --keep_non_rsid \
        --overwrite \
        --additional_output \
        --model_db_snp_key varID \
        --throw \
        --output_file /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/${tissue_name}/${ct_name}/${tissue_name}_${ct_name}_${trait_name}
    done
  done
done

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID