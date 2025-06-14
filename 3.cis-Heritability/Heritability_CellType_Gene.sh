#!/bin/bash

cd /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/Cell_Heritability/Cell_expression/MammaryGland/
phenotype_file="MammaryGland_Bulk.txt"
gene_list=$(head -n 1 $phenotype_file | awk '{for(i=3;i<=NF;i++) print $i}')
for gene in $gene_list; do
    gcta64 --bfile "gene_snp/${gene}_snps" --make-grm --out "grm/${gene}_grm"
    
    gcta64 --grm "grm/${gene}_grm" \
           --pheno "Bulk_output/${gene}_pheno.txt" \
           --qcovar MammaryGland_cov.txt \
           --reml-est-fix \
           --reml \
           --out "Bulk_output/${gene}_heritability"
done
