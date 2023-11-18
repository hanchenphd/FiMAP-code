##### randomly select ultra-rare causal variants and assign effect sizes
### single segment analyzed using ibd_mapping_batch
Sys.setenv(MKL_NUM_THREADS = 1)
library(GMMAT)
library(Matrix)
library(doMC)
library(FiMAP)
seed <- 1
set.seed(seed)
ids <- read.table("ids_after_filtering")[,1]
ids_exclude <- read.table("withdrawal.csv")[,1]
ids <- ids[!ids %in% ids_exclude]
N <- length(ids)
print(N)
plink.ids <- read.table("ukb_hap_merged.fam")[,2]
id.include <- which(ids %in% plink.ids)
N.randomvec <- 100
ncores <- 1
nsel <- 2
rsq <- 0.02
cat("Effect size rsq:", rsq, "\n")
maf1 <- 0.0005
maf2 <- 0.01
info <- read.table("geno_counts.csv", sep = ",")
info <- subset(info, V5 >= nsel) # 0.0005
ibd3 <- get(load("3cm/1/all_3.max.global.cut0.088.RData"))
ibd3 <- ibd3[rownames(ibd3) %in% ids, colnames(ibd3) %in% ids]
ibd5 <- get(load("5cm/1/all_5.max.global.cut0.088.RData"))
ibd5 <- ibd5[rownames(ibd5) %in% ids, colnames(ibd5) %in% ids]
ibd10 <- get(load("10cm/1/all_10.max.global.cut0.088.RData"))
ibd10 <- ibd10[rownames(ibd10) %in% ids, colnames(ibd10) %in% ids]
hapIBD.ibd3 <- get(load("hapIBD/3cm/all_3.ibd.global.cut0.088.RData"))
hapIBD.ibd3 <- hapIBD.ibd3[rownames(hapIBD.ibd3) %in% ids, colnames(hapIBD.ibd3) %in% ids]
hapIBD.ibd5 <- get(load("hapIBD/5cm/all_5.ibd.global.cut0.088.RData"))
hapIBD.ibd5 <- hapIBD.ibd5[rownames(hapIBD.ibd5) %in% ids, colnames(hapIBD.ibd5) %in% ids]
hapIBD.ibd5 <- hapIBD.ibd5 + 0.002 * Diagonal(n = nrow(hapIBD.ibd5)) # positive definite
hapIBD.ibd10 <- get(load("hapIBD/10cm/all_10.ibd.global.cut0.088.RData"))
hapIBD.ibd10 <- hapIBD.ibd10[rownames(hapIBD.ibd10) %in% ids, colnames(hapIBD.ibd10) %in% ids]
IBD.chol <- chol(ibd3) # used in the generative model

seed <- as.numeric(commandArgs(T)[1]) # 1-40
cat("Using seed", seed, "\n")
set.seed(seed)
nsimul <- 25
res <- NULL
for(i in 1:nsimul) {
    print(i)
    ii <- sample(nrow(info), 1)
    print(info[ii,])
    chr <- info[ii,1]
    start <- info[ii,2]
    start1 <- start - 1
    start2 <- start + 1
    ibd3.file <- paste0("3cm/1/chunk.1cM/chr",chr,".start",start,"cM.RData")
    ibd5.file <- paste0("5cm/1/chunk.1cM/chr",chr,".start",start,"cM.RData")
    ibd10.file <- paste0("10cm/1/chunk.1cM/chr",chr,".start",start,"cM.RData")
    hapIBD.ibd3.file <- paste0("hapIBD/3cm/chunk.1cM/chr",chr,".start",start,"cM.RData")
    hapIBD.ibd5.file <- paste0("hapIBD/5cm/chunk.1cM/chr",chr,".start",start,"cM.RData")
    hapIBD.ibd10.file <- paste0("hapIBD/10cm/chunk.1cM/chr",chr,".start",start,"cM.RData")
    gds.file <- paste0("data/UKB/plink/chunk_5.1cM/chr",chr,".start",start,"cM.gds")
    gds2.file <- paste0("data/UKB/imputed/cm/chunk_1cM/chr",chr,".start",start,"cM.gds")
    ibd3.file1 <- paste0("3cm/1/chunk.1cM/chr",chr,".start",start1,"cM.RData")
    ibd5.file1 <- paste0("5cm/1/chunk.1cM/chr",chr,".start",start1,"cM.RData")
    ibd10.file1 <- paste0("10cm/1/chunk.1cM/chr",chr,".start",start1,"cM.RData")
    hapIBD.ibd3.file1 <- paste0("hapIBD/3cm/chunk.1cM/chr",chr,".start",start1,"cM.RData")
    hapIBD.ibd5.file1 <- paste0("hapIBD/5cm/chunk.1cM/chr",chr,".start",start1,"cM.RData")
    hapIBD.ibd10.file1 <- paste0("hapIBD/10cm/chunk.1cM/chr",chr,".start",start1,"cM.RData")
    gds.file1 <- paste0("data/UKB/plink/chunk_5.1cM/chr",chr,".start",start1,"cM.gds")
    file1 <- file.exists(gds.file1)
    gds2.file1 <- paste0("data/UKB/imputed/cm/chunk_1cM/chr",chr,".start",start1,"cM.gds")
    ibd3.file2 <- paste0("3cm/1/chunk.1cM/chr",chr,".start",start2,"cM.RData")
    ibd5.file2 <- paste0("5cm/1/chunk.1cM/chr",chr,".start",start2,"cM.RData")
    ibd10.file2 <- paste0("10cm/1/chunk.1cM/chr",chr,".start",start2,"cM.RData")
    hapIBD.ibd3.file2 <- paste0("hapIBD/3cm/chunk.1cM/chr",chr,".start",start2,"cM.RData")
    hapIBD.ibd5.file2 <- paste0("hapIBD/5cm/chunk.1cM/chr",chr,".start",start2,"cM.RData")
    hapIBD.ibd10.file2 <- paste0("hapIBD/10cm/chunk.1cM/chr",chr,".start",start2,"cM.RData")
    gds.file2 <- paste0("data/UKB/plink/chunk_5.1cM/chr",chr,".start",start2,"cM.gds")
    file2 <- file.exists(gds.file2)
    gds2.file2 <- paste0("data/UKB/imputed/cm/chunk_1cM/chr",chr,".start",start2,"cM.gds")
    gds <- SeqArray::seqOpen(gds.file)
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    sample.id <- sample.id[sample.id %in% ids]
    SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
    geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
    SeqArray::seqClose(gds)
    freq <- colMeans(geno, na.rm = T)/2
    beta <- ifelse(freq>0 & freq<1 & (freq<maf1 | freq>1-maf1), sqrt(rsq/2/freq/(1-freq)), 0)
    candidates <- which(beta != 0)
    selected <- sort(sample(candidates, nsel))
    geneffects <- as.numeric(geno[,selected] %*% beta[selected])

    go_on <- TRUE
    while(go_on) {
        age <- rnorm(N, 50, 5)
	sex <- rbinom(N, 1, 0.5)
	raneff <- as.numeric(crossprod(IBD.chol, rnorm(N)))
	ranerr <- rnorm(N)
	y <- 0.05*age+0.5*sex+raneff+ranerr
	pheno <- data.frame(id = ids, y = y, age = age, sex = sex)
	pheno <- pheno[sample(id.include, 400000),]
	pheno$y <- pheno$y + geneffects[match(pheno$id, sample.id)]
	obj3 <- glmmkin(y ~ age + sex, data = pheno, id = "id", kins = ibd3, family = gaussian(link = "identity"), verbose = F)
	obj5 <- glmmkin(y ~ age + sex, data = pheno, id = "id", kins = ibd5, family = gaussian(link = "identity"), verbose = F)
	obj10 <- glmmkin(y ~ age + sex, data = pheno, id = "id", kins = ibd10, family = gaussian(link = "identity"), verbose = F)
	hapIBD.obj3 <- glmmkin(y ~ age + sex, data = pheno, id = "id", kins = hapIBD.ibd3, family = gaussian(link = "identity"), verbose = F)
	hapIBD.obj5 <- glmmkin(y ~ age + sex, data = pheno, id = "id", kins = hapIBD.ibd5, family = gaussian(link = "identity"), verbose = F)
	hapIBD.obj10 <- glmmkin(y ~ age + sex, data = pheno, id = "id", kins = hapIBD.ibd10, family = gaussian(link = "identity"), verbose = F)
        if(all(c(obj3$theta,obj5$theta,obj10$theta,hapIBD.obj3$theta,hapIBD.obj5$theta,hapIBD.obj10$theta)>0)) {
	    go_on <- FALSE
	}
    }
    match.idx1 <- match(pheno$id, rownames(ibd5))
    match.idx2 <- match(pheno$id, colnames(ibd5))
    tmpibd <- ibd5[match.idx1, match.idx2]
    tmpibd.chol <- chol(tmpibd)
    Sigma <- obj5$theta[1]*Diagonal(nrow(tmpibd))+obj5$theta[2]*tmpibd
    OBJ <- glmmkin2randomvec(obj5, Z = list(t(tmpibd.chol)), N.randomvec = N.randomvec)
    offset <- crossprod(OBJ$random.vectors, crossprod(Sigma, OBJ$random.vectors))
    out5 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd5.file)
    out5$n.nonzero5 <- out5$n.nonzero
    out5$p.value5 <- out5$p.value
    if(file1) {
    out5_1 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd5.file1)
    out5_1$n.nonzero5_1 <- out5_1$n.nonzero
    out5_1$p.value5_1 <- out5_1$p.value
    } else {
    out5_1 <- data.frame(n.nonzero5_1=NA, p.value5_1=NA)
    }
    if(file2) {
    out5_2 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd5.file2)
    out5_2$n.nonzero5_2 <- out5_2$n.nonzero
    out5_2$p.value5_2 <- out5_2$p.value
    } else {
    out5_2 <- data.frame(n.nonzero5_2=NA, p.value5_2=NA)
    }
    match.idx1 <- match(pheno$id, rownames(ibd3))
    match.idx2 <- match(pheno$id, colnames(ibd3))
    tmpibd <- ibd3[match.idx1, match.idx2]
    tmpibd.chol <- chol(tmpibd)
    Sigma <- obj3$theta[1]*Diagonal(nrow(tmpibd))+obj3$theta[2]*tmpibd
    OBJ <- glmmkin2randomvec(obj3, Z = list(t(tmpibd.chol)), N.randomvec = N.randomvec)
    offset <- crossprod(OBJ$random.vectors, crossprod(Sigma, OBJ$random.vectors))
    out3 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd3.file)
    out3$n.nonzero3 <- out3$n.nonzero
    out3$p.value3 <- out3$p.value
    if(file1) {
    out3_1 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd3.file1)
    out3_1$n.nonzero3_1 <- out3_1$n.nonzero
    out3_1$p.value3_1 <- out3_1$p.value
    } else {
    out3_1 <- data.frame(n.nonzero3_1=NA, p.value3_1=NA)
    }
    if(file2) {
    out3_2 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd3.file2)
    out3_2$n.nonzero3_2 <- out3_2$n.nonzero
    out3_2$p.value3_2 <- out3_2$p.value
    } else {
    out3_2 <- data.frame(n.nonzero3_2=NA, p.value3_2=NA)
    }
    match.idx1 <- match(pheno$id, rownames(ibd10))
    match.idx2 <- match(pheno$id, colnames(ibd10))
    tmpibd <- ibd10[match.idx1, match.idx2]
    tmpibd.chol <- chol(tmpibd)
    Sigma <- obj10$theta[1]*Diagonal(nrow(tmpibd))+obj10$theta[2]*tmpibd
    OBJ <- glmmkin2randomvec(obj10, Z = list(t(tmpibd.chol)), N.randomvec = N.randomvec)
    offset <- crossprod(OBJ$random.vectors, crossprod(Sigma, OBJ$random.vectors))
    out10 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd10.file)
    out10$n.nonzero10 <- out10$n.nonzero
    out10$p.value10 <- out10$p.value
    if(file1) {
    out10_1 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd10.file1)
    out10_1$n.nonzero10_1 <- out10_1$n.nonzero
    out10_1$p.value10_1 <- out10_1$p.value
    } else {
    out10_1 <- data.frame(n.nonzero10_1=NA, p.value10_1=NA)
    }
    if(file2) {
    out10_2 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = ibd10.file2)
    out10_2$n.nonzero10_2 <- out10_2$n.nonzero
    out10_2$p.value10_2 <- out10_2$p.value
    } else {
    out10_2 <- data.frame(n.nonzero10_2=NA, p.value10_2=NA)
    }
    match.idx1 <- match(pheno$id, rownames(hapIBD.ibd5))
    match.idx2 <- match(pheno$id, colnames(hapIBD.ibd5))
    tmpibd <- hapIBD.ibd5[match.idx1, match.idx2]
    tmpibd.chol <- chol(tmpibd)
    Sigma <- hapIBD.obj5$theta[1]*Diagonal(nrow(tmpibd))+hapIBD.obj5$theta[2]*tmpibd
    OBJ <- glmmkin2randomvec(hapIBD.obj5, Z = list(t(tmpibd.chol)), N.randomvec = N.randomvec)
    offset <- crossprod(OBJ$random.vectors, crossprod(Sigma, OBJ$random.vectors))
    hapIBD.out5 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd5.file)
    hapIBD.out5$n.nonzero.hapIBD5 <- hapIBD.out5$n.nonzero
    hapIBD.out5$p.value.hapIBD5 <- hapIBD.out5$p.value
    if(file1) {
    hapIBD.out5_1 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd5.file1)
    hapIBD.out5_1$n.nonzero.hapIBD5_1 <- hapIBD.out5_1$n.nonzero
    hapIBD.out5_1$p.value.hapIBD5_1 <- hapIBD.out5_1$p.value
    } else {
    hapIBD.out5_1 <- data.frame(n.nonzero.hapIBD5_1=NA, p.value.hapIBD5_1=NA)
    }
    if(file2) {
    hapIBD.out5_2 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd5.file2)
    hapIBD.out5_2$n.nonzero.hapIBD5_2 <- hapIBD.out5_2$n.nonzero
    hapIBD.out5_2$p.value.hapIBD5_2 <- hapIBD.out5_2$p.value
    } else {
    hapIBD.out5_2 <- data.frame(n.nonzero.hapIBD5_2=NA, p.value.hapIBD5_2=NA)
    }
    match.idx1 <- match(pheno$id, rownames(hapIBD.ibd3))
    match.idx2 <- match(pheno$id, colnames(hapIBD.ibd3))
    tmpibd <- hapIBD.ibd3[match.idx1, match.idx2]
    tmpibd.chol <- chol(tmpibd)
    Sigma <- hapIBD.obj3$theta[1]*Diagonal(nrow(tmpibd))+hapIBD.obj3$theta[2]*tmpibd
    OBJ <- glmmkin2randomvec(hapIBD.obj3, Z = list(t(tmpibd.chol)), N.randomvec = N.randomvec)
    offset <- crossprod(OBJ$random.vectors, crossprod(Sigma, OBJ$random.vectors))
    hapIBD.out3 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd3.file)
    hapIBD.out3$n.nonzero.hapIBD3 <- hapIBD.out3$n.nonzero
    hapIBD.out3$p.value.hapIBD3 <- hapIBD.out3$p.value
    if(file1) {
    hapIBD.out3_1 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd3.file1)
    hapIBD.out3_1$n.nonzero.hapIBD3_1 <- hapIBD.out3_1$n.nonzero
    hapIBD.out3_1$p.value.hapIBD3_1 <- hapIBD.out3_1$p.value
    } else {
    hapIBD.out3_1 <- data.frame(n.nonzero.hapIBD3_1=NA, p.value.hapIBD3_1=NA)
    }
    if(file2) {
    hapIBD.out3_2 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd3.file2)
    hapIBD.out3_2$n.nonzero.hapIBD3_2 <- hapIBD.out3_2$n.nonzero
    hapIBD.out3_2$p.value.hapIBD3_2 <- hapIBD.out3_2$p.value
    } else {
    hapIBD.out3_2 <- data.frame(n.nonzero.hapIBD3_2=NA, p.value.hapIBD3_2=NA)
    }
    match.idx1 <- match(pheno$id, rownames(hapIBD.ibd10))
    match.idx2 <- match(pheno$id, colnames(hapIBD.ibd10))
    tmpibd <- hapIBD.ibd10[match.idx1, match.idx2]
    tmpibd.chol <- chol(tmpibd)
    Sigma <- hapIBD.obj10$theta[1]*Diagonal(nrow(tmpibd))+hapIBD.obj10$theta[2]*tmpibd
    OBJ <- glmmkin2randomvec(hapIBD.obj10, Z = list(t(tmpibd.chol)), N.randomvec = N.randomvec)
    offset <- crossprod(OBJ$random.vectors, crossprod(Sigma, OBJ$random.vectors))
    hapIBD.out10 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd10.file)
    hapIBD.out10$n.nonzero.hapIBD10 <- hapIBD.out10$n.nonzero
    hapIBD.out10$p.value.hapIBD10 <- hapIBD.out10$p.value
    if(file1) {
    hapIBD.out10_1 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd10.file1)
    hapIBD.out10_1$n.nonzero.hapIBD10_1 <- hapIBD.out10_1$n.nonzero
    hapIBD.out10_1$p.value.hapIBD10_1 <- hapIBD.out10_1$p.value
    } else {
    hapIBD.out10_1 <- data.frame(n.nonzero.hapIBD10_1=NA, p.value.hapIBD10_1=NA)
    }
    if(file2) {
    hapIBD.out10_2 <- ibd_mapping(OBJ, offset = offset, ibd.ids = ids, ibd.file = hapIBD.ibd10.file2)
    hapIBD.out10_2$n.nonzero.hapIBD10_2 <- hapIBD.out10_2$n.nonzero
    hapIBD.out10_2$p.value.hapIBD10_2 <- hapIBD.out10_2$p.value
    } else {
    hapIBD.out10_2 <- data.frame(n.nonzero.hapIBD10_2=NA, p.value.hapIBD10_2=NA)
    }
    out <- cbind(out3[c("chr", "start", "end", "n.nonzero3", "p.value3")], out3_1[c("n.nonzero3_1", "p.value3_1")], out3_2[c("n.nonzero3_2", "p.value3_2")], out5[c("n.nonzero5", "p.value5")], out5_1[c("n.nonzero5_1", "p.value5_1")], out5_2[c("n.nonzero5_2", "p.value5_2")], out10[c("n.nonzero10", "p.value10")], out10_1[c("n.nonzero10_1", "p.value10_1")], out10_2[c("n.nonzero10_2", "p.value10_2")], hapIBD.out3[c("n.nonzero.hapIBD3", "p.value.hapIBD3")], hapIBD.out3_1[c("n.nonzero.hapIBD3_1", "p.value.hapIBD3_1")], hapIBD.out3_2[c("n.nonzero.hapIBD3_2", "p.value.hapIBD3_2")], hapIBD.out5[c("n.nonzero.hapIBD5", "p.value.hapIBD5")], hapIBD.out5_1[c("n.nonzero.hapIBD5_1", "p.value.hapIBD5_1")], hapIBD.out5_2[c("n.nonzero.hapIBD5_2", "p.value.hapIBD5_2")], hapIBD.out10[c("n.nonzero.hapIBD10", "p.value.hapIBD10")], hapIBD.out10_1[c("n.nonzero.hapIBD10_1", "p.value.hapIBD10_1")], hapIBD.out10_2[c("n.nonzero.hapIBD10_2", "p.value.hapIBD10_2")])
    glmm.score(obj3, infile = gds.file, outfile = paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start,"cM.txt"), MAF.range = c(maf2, 0.5), miss.cutoff = 0.05, ncores = ncores)
    svtres <- read.table(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start,"cM.txt"), header=T, as.is=T)
    out$p.svt.min <- min(svtres$PVAL)
    out$p.svt.cauchy <- pcauchy(mean(qcauchy(svtres$PVAL,lower=F)),lower=F)
    unlink(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start,"cM.txt"))
    glmm.score(obj3, infile = gds2.file, outfile = paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start,"cM.txt"), MAF.range = c(maf2, 0.5), miss.cutoff = 0.05, ncores = ncores)
    svtres <- read.table(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start,"cM.txt"), header=T, as.is=T)
    out$p.svt.imputed.min <- min(svtres$PVAL)
    out$p.svt.imputed.cauchy <- pcauchy(mean(qcauchy(svtres$PVAL,lower=F)),lower=F)
    unlink(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start,"cM.txt"))
    if(file1) {
    glmm.score(obj3, infile = gds.file1, outfile = paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start1,"cM.txt"), MAF.range = c(maf2, 0.5), miss.cutoff = 0.05, ncores = ncores)
    svtres <- read.table(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start1,"cM.txt"), header=T, as.is=T)
    out$p.svt.min1 <- min(svtres$PVAL)
    out$p.svt.cauchy1 <- pcauchy(mean(qcauchy(svtres$PVAL,lower=F)),lower=F)
    unlink(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start1,"cM.txt"))
    glmm.score(obj3, infile = gds2.file1, outfile = paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start1,"cM.txt"), MAF.range = c(maf2, 0.5), miss.cutoff = 0.05, ncores = ncores)
    svtres <- read.table(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start1,"cM.txt"), header=T, as.is=T)
    out$p.svt.imputed.min1 <- min(svtres$PVAL)
    out$p.svt.imputed.cauchy1 <- pcauchy(mean(qcauchy(svtres$PVAL,lower=F)),lower=F)
    unlink(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start1,"cM.txt"))
    } else out$p.svt.min1 <- out$p.svt.cauchy1 <- out$p.svt.imputed.min1 <- out$p.svt.imputed.cauchy1 <- NA
    if(file2) {
    glmm.score(obj3, infile = gds.file2, outfile = paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start2,"cM.txt"), MAF.range = c(maf2, 0.5), miss.cutoff = 0.05, ncores = ncores)
    svtres <- read.table(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start2,"cM.txt"), header=T, as.is=T)
    out$p.svt.min2 <- min(svtres$PVAL)
    out$p.svt.cauchy2 <- pcauchy(mean(qcauchy(svtres$PVAL,lower=F)),lower=F)
    unlink(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start2,"cM.txt"))
    glmm.score(obj3, infile = gds2.file2, outfile = paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start2,"cM.txt"), MAF.range = c(maf2, 0.5), miss.cutoff = 0.05, ncores = ncores)
    svtres <- read.table(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start2,"cM.txt"), header=T, as.is=T)
    out$p.svt.imputed.min2 <- min(svtres$PVAL)
    out$p.svt.imputed.cauchy2 <- pcauchy(mean(qcauchy(svtres$PVAL,lower=F)),lower=F)
    unlink(paste0("untyped/rsq_",rsq,"_seed", seed, ".GWAS.simul",i,".chr",chr,".start",start2,"cM.txt"))
    } else out$p.svt.min2 <- out$p.svt.cauchy2 <- out$p.svt.imputed.min2 <- out$p.svt.imputed.cauchy2 <- NA
    res <- rbind(res, out)
}
write.table(res, paste0("untyped/rsq_",rsq,"_seed", seed, ".txt"), quote=F, row.names=F, col.names=T)
