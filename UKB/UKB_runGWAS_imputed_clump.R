dict <- cbind(c(2:4, 35:37), c(48:50, 20015, 21001, 21002))
cmd <- commandArgs(T)
pheno.idx <- as.numeric(cmd[1]) # 48, 49, 50, 20015, 21001, 21002
cat("Analyzing phenotype", pheno.idx, "...\n")
maf <- as.numeric(cmd[2]) # 0.0001
cat("Use MAF filter", maf, "...\n")
flank <- as.numeric(cmd[3]) # 100
cat("Clumping with flanking region", flank, "kb...\n")
alpha <- as.numeric(cmd[4]) # 1e-6
cat("GWAS SNP p-value threshold at", alpha, "...\n")
qual <- 0.3 # imputation quality filter
missrate <- 0.05

for(chr in 1:22) {    
  gwas <- read.table(paste0("GWAS_imputed.chr",chr,".X", pheno.idx, ".0.0.QUALgt0.3.txt"), header = T, as.is = T)
  gwas <- subset(gwas, AF >= maf & AF <= 1 - maf)
  gwas <- subset(gwas, MISSRATE < missrate)
  gwas <- subset(gwas, QUAL > qual) # imputed
  gwas <- subset(gwas, !is.na(PVAL) & PVAL < alpha)
  curr_pval <- min(gwas$PVAL)
  while(T) {
    curr_pval_idx <- which(gwas$PVAL==curr_pval)
    idx <- c()
    for(curr_idx in curr_pval_idx) {
      if(curr_idx %in% idx) next
      idx <- unique(c(idx, setdiff(which(gwas$PVAL>=curr_pval & abs(gwas$POS-gwas$POS[curr_idx])<=flank*1000), curr_idx)))
    }
    if(length(idx)>0) gwas <- gwas[-idx,]
    if(all(gwas$PVAL<=curr_pval)) break
    curr_pval <- min(gwas$PVAL[gwas$PVAL>curr_pval])
  }
  outfile <- paste0("GWAS_imputed.chr",chr,".pheno", pheno.idx, "_maf",maf,"_missrate",missrate,"_qual",qual,"_clump",flank,"kb_alpha",alpha,".out")
  write.table(gwas, outfile, quote = F, row.names = F, col.names = T, sep = "\t")
}
