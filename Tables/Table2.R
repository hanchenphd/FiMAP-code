# Table 2. Numbers of significant 1 cM windows in UK Biobank unconditional and conditional IBD mapping for 6 anthropometric traits.

# Waist circumference (48.0.0)
# Hip circumference (49.0.0)
# Standing height (50.0.0)
# Sitting height (20015.0.0)
# Body mass index (21001.0.0)
# Weight (21002.0.0)

library(data.table)
maf <- 0.0001
missrate <- 0.05
qual <- 0.3
pgwas <- 1e-6
pfimap <- 0.05/3403
ct <- ct1 <- ct2 <- 0
for(pheno in c(48,49,50,20015,21001,21002)) {
data1 <- NULL
for(CHR in 1:22) {
infile <- paste0("GWAS_imputed.chr",CHR,".X",pheno,".0.0.QUALgt0.3.txt")
data <- fread(infile, data.table = F)
data <- subset(data, AF >= maf & AF <= 1 - maf)
data <- subset(data, MISSRATE < missrate)
data <- subset(data, QUAL > qual) # imputed
data <- subset(data, !is.na(PVAL) & PVAL < pgwas)
data$start <- floor(data$cM)
data <- data[,c("CHR", "start", "PVAL")]
data1 <- rbind(data1, data)
}
# RaPID
infile <- paste0("3cm/1/seed12345_X",pheno,".0.0.txt")
data2 <- read.table(infile, header=T)
data2 <- data2[, c("chr", "start", "p.value")]
data2 <- subset(data2, p.value < pfimap)
data1$region <- paste(data1$CHR, data1$start, sep = ":")
data2$region <- paste(data2$chr, data2$start, sep = ":")
cat("Phenotype", pheno, "...\n")
cat("RaPID FiMAP:", nrow(data2), "hits...\n")
cat("No GWAS overlaps:", sum(!data2$region %in% data1$region), "hits...\n")
cat("Proportion:", sum(!data2$region %in% data1$region)/nrow(data2), "\n")
ct <- ct + nrow(data2)
ct1 <- ct1 + sum(!data2$region %in% data1$region)
infile <- paste0("3cm/1/conditional_imputed_GDS_pheno",pheno,"_maf1e-04_missrate0.05_qual0.3_flank3cM_alpha1e-06.20230714.out")
data3 <- read.table(infile, header=T)
data3 <- subset(data3, p.value.conditional < pfimap)
cat("RaPID conditional FiMAP:", nrow(data3), "hits...\n")
cat("Proportion:", nrow(data3)/nrow(data2), "\n")
ct2 <- ct2 + nrow(data3)
}
cat("\nSummary:\n")
cat("RaPID FiMAP:", ct, "total hits...\n")
cat("No GWAS overlaps:", ct1, "total hits...\n")
cat("Proportion:", ct1/ct, "\n")
cat("RaPID FiMAP conditional:", ct2, "total hits...\n")
cat("Proportion:", ct2/ct, "\n")
