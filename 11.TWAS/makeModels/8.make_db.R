library("RSQLite")
library("dplyr")

argv <- commandArgs(trailingOnly = TRUE)
type = "eQTL"
tissue <- argv[1]
celltype <- argv[2]

"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')

model_summaries <- read.table(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/result/chr1/eQTL.", tissue, "_", celltype, "_Model_training_chr1_model_summaries.txt"),header = T, stringsAsFactors = F)
tiss_summary <- read.table(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/result/chr1/eQTL.", tissue, "_", celltype, "_Model_training_chr1_summary.txt"), header = T, stringsAsFactors = F)
n_samples <- tiss_summary$n_samples
  
for (i in c(2:29)) {
  if (file.exists(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/result/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_model_summaries.txt"))) {
    if (!file.info(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/result/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_model_summaries.txt"))$size == 0) {
      if (!file.info(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/weights/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_weights.txt"))$size == 0) {
        model_summaries <- rbind(model_summaries,
                            read.table(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/result/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_model_summaries.txt"), header = T, stringsAsFactors = F))
        tiss_summary <- rbind(tiss_summary,
                             read.table(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/result/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_summary.txt"), header = T, stringsAsFactors = F))
      }
    }                            
  }
}
  
model_summaries <- rename(model_summaries, gene = gene_id)

# Create a database connection
if (!dir.exists(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/dbs/"))) {
  dir.create(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/dbs/"))
}
conn <- dbConnect(drv = driver, paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/dbs/eQTL_", tissue, "_", celltype, "_ElasticNet_models.db"))
dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")

# Weights Table -----
weights <- read.table(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/weights/chr1/eQTL.", tissue, "_", celltype, "_Model_training_chr1_weights.txt"), header = T,stringsAsFactors = F)

for (i in c(2:29)) {
  if (file.exists(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/weights/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_weights.txt"))) {
    if (!file.info(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/weights/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_weights.txt"))$size == 0) {
    weights <- rbind(weights,
                            read.table(paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/weights/chr", i, "/eQTL.", tissue, "_", celltype, "_Model_training_chr", i, "_weights.txt"), header = T, stringsAsFactors = F))
    }
  }
}
  
weights <- rename(weights, gene = gene_id)
weights$varID <- weights$rsid
dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

# Sample_info Table ----
sample_info <- data.frame(n_samples = n_samples, population = 'EUR') # Provide the population info
dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
  
# Construction Table ----
construction <- tiss_summary %>%
                    select(chrom, cv_seed) %>%
                    rename(chromosome = chrom)

dbWriteTable(conn, 'construction', construction, overwrite = TRUE)
dbDisconnect(conn)

#####
unfiltered_db <- paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/dbs/eQTL_", tissue, "_", celltype, "_ElasticNet_models.db")
filtered_db <- paste0("/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/TWAS/ct_seQTL/", tissue, "/", celltype, "/dbs/eQTL_", tissue, "_", celltype, "_ElasticNet_models_filter_0.05.db")
driver <- dbDriver("SQLite")
in_conn <- dbConnect(driver, unfiltered_db)
out_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from model_summaries where zscore_pval < 0.05 and rho_avg > 0.05')
model_summaries <- model_summaries %>% 
                    rename(pred.perf.R2 = rho_avg_squared, genename = gene_name, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model)
if (nrow(model_summaries)>0) {
  model_summaries$pred.perf.qval <- NA
  dbWriteTable(out_conn, 'extra', model_summaries, overwrite = TRUE)
  construction <- dbGetQuery(in_conn, 'select * from construction')
  dbWriteTable(out_conn, 'construction', construction, overwrite = TRUE)
  sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
  dbWriteTable(out_conn, 'sample_info', sample_info, overwrite = TRUE)
  weights <- dbGetQuery(in_conn, 'select * from weights')
  weights <- weights %>% filter(gene %in% model_summaries$gene) %>% rename(eff_allele = alt, ref_allele = ref, weight = beta)
  dbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)
  dbExecute(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbExecute(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbDisconnect(in_conn)
  dbDisconnect(out_conn)
}

