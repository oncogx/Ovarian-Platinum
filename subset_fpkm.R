setwd('/Volumes/oncogxA/Projects/PCAWG')
RNA_SEQ_SUBSET_PATH = "/Volumes/oncogxA/Projects/PCAWG/RNA/match_rna_by_icgc.tsv"

specimens <- read.csv(RNA_SEQ_SUBSET_PATH,sep="\t", header = TRUE)
ids <- specimens$aliquot_id

fpkm_full <- read.csv("RNA/joint_fpkm.tsv",sep="\t", header = TRUE)
fpkm_subset = fpkm_full[,ids]

write.csv(fpkm_subset, "RNA/joint_fpkm_subset_rna.csv")