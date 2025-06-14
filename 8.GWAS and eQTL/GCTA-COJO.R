path <- "/faststorage/project/cattle_gtexs/Downstream_analysis/complextraits/GWAS/" ##change to your GWAS summary data path
trait <- list.dirs(path, full.names = TRUE, recursive = FALSE)
trait <- sapply(trait, function(x) unlist(strsplit(x, "\\/"))[9])
trait <- data.frame(trait)
trait <- trait$trait

for (i in 1:length(trait)) {
  gwas <- fread(paste0(path, trait[[i]], "/", trait[[i]], "_china.mlma"))
  gwas <- data.frame(gwas)
  gwas$rs_id <- paste0(gwas[,2], "_", gwas[,3])
  select_gwas <- gwas[,c(13,6,5,7,8,9,10,12)]
  names(select_gwas) <- c("SNP","A1","A2","freq","b","se","p","N")
  write.table(select_gwas, paste0(path, trait[[i]], "/", trait[[i]], "_cojo.txt"), sep = "\t", row.names = F, quote = FALSE)
}

main_dir="/faststorage/project/cattle_gtexs/Downstream_analysis/complextraits/GWAS/" ##change to your GWAS summary data path
bfile_file="/faststorage/project/cattle_gtexs/CattleGTEx/panel_Hols" ##change to your reference panel path
for trait_dir in "$main_dir"/*/; do
  gwas_file="${trait_dir}$(basename "$trait_dir")_cojo.txt"
  out_dir="${trait_dir}$(basename "$trait_dir")_cojo_res"
  if [ -e "$gwas_file" ]; then
  gcta64 --bfile "$bfile_file" \
         --maf 0.05 \
         --cojo-file "$gwas_file" \
         --cojo-slct \
         --cojo-p 1e-5 \
         --out "$out_dir"
  fi
done