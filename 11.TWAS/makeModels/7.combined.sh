#!/bin/bash
work_path=$1
type=$2
tissue=$3

 
 
base_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
for tissue_dir in "$base_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  for tt in {1..29}; do
    cd /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/bulk_eQTL/${tissue_name}/result/chr${tt}/
    cat eQTL.${tissue_name}_Model_training_chr${tt}.*.model_summaries.txt > eQTL.${tissue_name}_Model_training_chr${tt}_model_summaries.txt
    sed -i '1i\gene_id\tgene_name\tgene_type\talpha\tn_snps_in_window\tn_snps_in_model\tlambda_min_mse\ttest_R2_avg\ttest_R2_sd\tcv_R2_avg\tcv_R2_sd\tin_sample_R2\tnested_cv_fisher_pval\trho_avg\trho_se\trho_zscore\trho_avg_squared\tzscore_pval\tcv_rho_avg\tcv_rho_se\tcv_rho_avg_squared\tcv_zscore_est\tcv_zscore_pval\tcv_pval_est' eQTL.${tissue_name}_Model_training_chr${tt}_model_summaries.txt
  done
done

base_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
for tissue_dir in "$base_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  for ttt in {1..29}; do
    cd /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/bulk_eQTL/${tissue_name}/covariances/chr${ttt}/
    cat eQTL.${tissue_name}_Model_training_chr${ttt}.*.covariances.txt > eQTL.${tissue_name}_Model_training_chr${ttt}_covariances.txt
    sed -i '1i\GENE RSID1 RSID2 VALUE' eQTL.${tissue_name}_Model_training_chr${ttt}_covariances.txt
  done
done

base_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
for tissue_dir in "$base_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  for tttt in {1..29}; do
    cd /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/bulk_eQTL/${tissue_name}/weights/chr${tttt}/
    cat eQTL.${tissue_name}_Model_training_chr${tttt}.*.weights.txt > eQTL.${tissue_name}_Model_training_chr${tttt}_weights.txt
    sed -i '1i\gene_id\trsid\tvarID\tref\talt\tbeta' eQTL.${tissue_name}_Model_training_chr${tttt}_weights.txt
  done
done
