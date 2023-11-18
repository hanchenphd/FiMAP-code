cmd <- commandArgs(T)
pheno.idx <- as.numeric(cmd[1]) # 48, 49, 50, 20015, 21001, 21002
chr <- as.numeric(cmd[2]) # 1-22
seed <- 12345
runs <- 1
qual <- as.numeric(cmd[3]) # 0.3

infile <- paste0("GWAS_imputed.chr",chr,".X",pheno.idx,".0.0.txt")
outfile <- paste0("GWAS_imputed.chr",chr,".X",pheno.idx,".0.0.QUALgt",qual,".txt")
infofile <- paste0("ukb_mfi_chr",chr,"_v3.txt")
mapfile <- paste0("ukb_",chr,".rMap")
bimfile <- "ukb_hap_merged.bim"
cat("Processing", infile, "...\n")

cm <- read.table(mapfile)[,2]
bim <- read.table(bimfile, as.is = T)
bim <- subset(bim, V1 == chr)
bim <- cbind(bim, cm)
info <- read.table(infofile,as.is=T)
info$SNP1 <- paste(chr, info$V3, info$V5, info$V4, sep=":")
tmp <- read.table(infile,header=T,as.is=T)
SNP1 <- paste(tmp$CHR, tmp$POS, tmp$REF, tmp$ALT, sep=":")
tmp$QUAL <- info$V8[match(SNP1, info$SNP1)]
tmp <- subset(tmp, QUAL > qual)
tmp$cM <- NA

i <- 1
for(j in 1:nrow(tmp)) {
    if(tmp$POS[j] < bim$V4[1]) {
        tmp$cM[j] <- bim$cm[1]
	next
    }
    while(tmp$POS[j] >= bim$V4[i] && i < nrow(bim)) i <- i + 1
    if(tmp$POS[j] < bim$V4[i]) tmp$cM[j] <- bim$cm[i-1] + (tmp$POS[j]-bim$V4[i-1])*(bim$cm[i]-bim$cm[i-1])/(bim$V4[i]-bim$V4[i-1])
    else {
        tmp$cM[j:nrow(tmp)] <- bim$cm[nrow(bim)]
	break
    }
}
write.table(tmp, outfile, quote = F, col.names = T, row.names = F, sep = "\t", na = "NA")
