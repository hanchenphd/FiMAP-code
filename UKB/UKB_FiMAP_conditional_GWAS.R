Sys.setenv(MKL_NUM_THREADS = 1)
library(GMMAT)
library(Matrix)
library(doMC)
library(FiMAP)

dict <- cbind(c(2:4, 35:37), c(48:50, 20015, 21001, 21002))
cmd <- commandArgs(T)
pheno.idx <- as.numeric(cmd[1]) # 48, 49, 50, 20015, 21001, 21002
cat("Analyzing phenotype", pheno.idx, "...\n")
maf <- 0.0001
cat("Use MAF filter", maf, "...\n")
flank <- 3
cat("Include flanking region", flank, "cM...\n")
alpha <- 1e-6
cat("GWAS SNP p-value threshold at", alpha, "...\n")
seed <- 12345
cm <- 3
cat("RaPID", cm, "cm...\n")
runs <- 1
cat("FiMap: seed", seed, "and run", runs, "...\n")
ncores <- as.numeric(cmd[2]) # 40
qual <- 0.3 # imputation quality filter
missrate <- 0.05

global.ibd.dir <- paste0(cm,"cm/",runs,"/")
local.ibd.dir <- paste0(cm,"cm/",runs,"/chunk.1cM/")
outfile.dir <- paste0(cm,"cm/",runs,"/")
geno.dir <- "data/UKB/imputed/"
fimap.dir <- paste0(cm,"cm/",runs,"/")

fimap <- read.table(paste0(fimap.dir, "seed",seed,"_X", pheno.idx, ".0.0.txt"), header = T, as.is = T)

cat("Use FiMap Bonferroni threshold at", 0.05/nrow(fimap), "...\n")
idx <- which(fimap$p.value < 0.05/nrow(fimap))
cat("Out of", nrow(fimap), "regions, found", length(idx), "candidate regions for conditional analysis...\n")
ncores <- min(ncores, length(idx), parallel::detectCores(logical = TRUE))
cat("Running conditional analysis using", ncores, "cores...\n")

### for conditional association analysis
rank_norm <- function(x) {
    naidx <- is.na(x)
    xx <- x[!naidx]
    n <- length(xx)
    x[!naidx] <- qnorm((rank(xx) - 0.5)/n)
    return(x)
}

pheno.file <- "data/UKB/pheno/ukb.baseline.WhiteBritish.qced.csv"
N.randomvec <- 100

ids <- read.table("ids_after_filtering")[,1]
ids_exclude <- read.table("withdrawal.csv")[,1]
ids <- ids[!ids %in% ids_exclude]
N <- length(ids)
print(N)

ibd <- get(load(paste0(global.ibd.dir, "all_",cm,".max.global.cut0.088.RData")))

out <- fimap[idx,]
names(out)[4:5] <- c("n.nonzero.unconditional", "p.value.unconditional")
out$SNPs <- NA
out$N <- NA
out$n.nonzero.conditional <- NA
out$p.value.conditional <- NA

pheno <- read.csv(pipe(paste0("cut -d, -f1,3-7,10-49,", 49+dict[dict[,2]==pheno.idx, 1], " ", pheno.file)), as.is = T)
pheno <- subset(pheno, eid %in% ids)
pheno <- pheno[!apply(is.na(pheno), 1, any), ]
pheno.name <- names(pheno)[length(names(pheno))]
cat("Phenotype colname:", pheno.name, "\n")
covars <- c("sex", "age", "age2", "age_sex", "age2_sex", paste0("PC", 1:10))
pheno$res <- rank_norm(residuals(lm(as.formula(paste(pheno.name, paste(covars, collapse = " + "), sep = " ~ ")), data = pheno)))
cat("Number of non-missing residuals:", sum(!is.na(pheno$res)), "\n")

doMC::registerDoMC(cores = ncores)
n.groups.percore <- (length(idx)-1) %/% ncores + 1
n.groups.percore_1 <- n.groups.percore * ncores - length(idx)
b <- NULL
tmp.out <- foreach(b=1:ncores, .combine=rbind, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
  tmp.idx <- if(b <= n.groups.percore_1) ((b-1)*(n.groups.percore-1)+1):(b*(n.groups.percore-1)) else (n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1-1)*n.groups.percore+1):(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1)*n.groups.percore)
  tmp2.out <- out[tmp.idx, , drop = FALSE]
  for(ii in tmp.idx) {
    cat("\n Analyze FiMap candidate region", ii, "...\n")
    ibd.file <- paste0(local.ibd.dir, "chr",out$chr[ii],".start",out$start[ii],"cM.RData")
    snp.selected <- c()
    
    gwas <- read.table(paste0("GWAS_imputed.chr",out$chr[ii],".pheno", pheno.idx, "_maf",maf,"_missrate",missrate,"_qual",qual,"_clump100kb_alpha",alpha,".out"), header = T, as.is = T)
    gwas$SNP <- paste(gwas$CHR, gwas$POS, gwas$REF, gwas$ALT, sep = ":")
    
    geno.file <- paste0(geno.dir, "chr",out$chr[ii],"_v3.impute2.gds")
    gds <- SeqArray::seqOpen(geno.file)
    GWASID <- SeqArray::seqGetData(gds, "sample.id")
    variant.id <- SeqArray::seqGetData(gds, "variant.id") # 1 2 3 4 5
    tmp.chr <- SeqArray::seqGetData(gds, "chromosome")
    tmp.pos <- SeqArray::seqGetData(gds, "position")
    alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
    tmp.ref <- unlist(lapply(alleles.list, function(x) x[1]))
    tmp.alt <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
    rm(alleles.list)
    snp.id <- paste(tmp.chr, tmp.pos, tmp.ref, tmp.alt, sep = ":")
    rm(tmp.chr, tmp.pos, tmp.ref, tmp.alt)

    tmp.pheno <- subset(pheno, eid %in% GWASID)
    GWASID <- GWASID[GWASID %in% tmp.pheno$eid]
    tmp.pheno.old <- tmp.pheno

    if(!any(gwas$CHR == out$chr[ii] & gwas$cM >= out$start[ii] - flank & gwas$cM < out$end[ii] + flank)) {
      SeqArray::seqClose(gds)
    } else {
      SeqArray::seqSetFilter(gds, sample.id = GWASID, variant.id = variant.id[snp.id %in% gwas$SNP[gwas$CHR == out$chr[ii] & gwas$cM >= out$start[ii] - flank & gwas$cM < out$end[ii] + flank]], verbose = TRUE)
      tmp.chr <- SeqArray::seqGetData(gds, "chromosome")
      tmp.pos <- SeqArray::seqGetData(gds, "position")
      alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
      tmp.ref <- unlist(lapply(alleles.list, function(x) x[1]))
      tmp.alt <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
      rm(alleles.list)
      snp.include <- paste(tmp.chr, tmp.pos, tmp.ref, tmp.alt, sep = ":")
      rm(tmp.chr, tmp.pos, tmp.ref, tmp.alt)

      snp.include2 <- gsub(",", ".", snp.include) # for R syntax compatibility
      snp.include2 <- gsub(":", ".", paste0("chr", snp.include2)) # for R syntax compatibility
      snp.selected <- snp.include
      snp.selected2 <- gsub(",", ".", snp.selected) # for R syntax compatibility
      snp.selected2 <- gsub(":", ".", paste0("chr", snp.selected2)) # for R syntax compatibility

      geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
      SeqArray::seqClose(gds)
      miss.idx <- which(is.na(geno))
      cat("Region", ii, ": Number of missing genotypes:", length(miss.idx), "\n")
      freq <- colMeans(geno, na.rm = TRUE)/2
      geno[miss.idx] <- 2*freq[ceiling(miss.idx/nrow(geno))]
      colnames(geno) <- snp.selected2
      tmp.pheno <- cbind(tmp.pheno.old, geno[match(tmp.pheno.old$eid, GWASID),,drop=FALSE])
    }
    if(length(snp.selected) > 0) {
      mdl <- lm(as.formula(paste("res ~", paste(covars, collapse = " + "), "+", paste(snp.selected2, collapse = " + "))), data = tmp.pheno)
      snp.pvals <- mdl$coef[(length(mdl$coef)-length(snp.selected)+1):length(mdl$coef)]
      if(any(is.na(snp.pvals))) {
        cat("Region", ii, ": Removing", snp.selected[which(is.na(snp.pvals))], "with NA LM p-value\n")
        snp.selected <- snp.selected[-which(is.na(snp.pvals))]
        snp.selected2 <- snp.selected2[-which(is.na(snp.pvals))]
      }
      tmp2.out$SNPs[match(ii, tmp.idx)] <- paste(snp.selected, collapse = ";")
      cat("Region", ii, ": Conditioning on SNPs", snp.selected, "\n")
      obj <- try(glmmkin(as.formula(paste("res ~", paste(snp.selected2, collapse = " + "), "+", paste(covars, collapse = " + "))), data = tmp.pheno, id = "eid", kins = ibd, family = gaussian(link = "identity")))
    } else {
      cat("Region", ii, ": Conditioning on NO SNPs\n")
      obj <- try(glmmkin(as.formula(paste("res ~", paste(covars, collapse = " + "))), data = tmp.pheno, id = "eid", kins = ibd, family = gaussian(link = "identity")))
    }
    if(class(obj) == "try-error") next
    cat("Region", ii, ": Final model: theta =", obj$theta, "\n")
    if(obj$theta[1]==0) next
    tmp2.out$N[match(ii, tmp.idx)] <- length(obj$id_include)
    set.seed(seed)
    match.idx1 <- match(pheno$eid, rownames(ibd))
    match.idx2 <- match(pheno$eid, colnames(ibd))
    tmpibd <- ibd[match.idx1, match.idx2]
    IBD.chol <- chol(tmpibd)
    Sigma <- obj$theta[1]*Diagonal(nrow(tmpibd))+obj$theta[2]*tmpibd
    obj1 <- glmmkin2randomvec(obj, Z = list(t(IBD.chol)), N.randomvec = N.randomvec, robust = FALSE)
    offset <- crossprod(obj1$random.vectors, crossprod(Sigma, obj1$random.vectors))
    rm(obj, tmpibd, Sigma)
    ibd.out <- ibd_mapping(obj1, offset = offset, ibd.ids = ids0, ibd.file = ibd.file)
    tmp2.out$n.nonzero.conditional[match(ii, tmp.idx)] <- ibd.out$n.nonzero
    tmp2.out$p.value.conditional[match(ii, tmp.idx)] <- ibd.out$p.value
  }
  tmp2.out
}
outfile <- paste0(outfile.dir, "conditional_imputed_GDS_pheno", pheno.idx, "_maf",maf,"_missrate",missrate,"_qual",qual,"_flank",flank,"cM_alpha",alpha,".out")
write.table(tmp.out, outfile, quote = F, row.names = F, col.names = T, sep = "\t")
