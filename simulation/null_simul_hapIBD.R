### single segment analyzed using ibd_mapping_batch
Sys.setenv(MKL_NUM_THREADS = 1)
library(GMMAT)
library(Matrix)
library(doMC)
library(FiMAP)
seed <- as.numeric(commandArgs(T)[1]) # 1-50
cm <- commandArgs(T)[2] # 3,5,10
runs <- commandArgs(T)[3] # 1
cat("Running null simulation using seed", seed, ", ", cm, "cm run", runs, "\n")
set.seed(seed)
ids <- read.table("ids_after_filtering")[,1]
ids_exclude <- read.table("withdrawal.csv")[,1]
ids <- ids[!ids %in% ids_exclude]
N <- length(ids)
print(N)
N.randomvec <- 100
ncores <- 40
ibd <- get(load(paste0("hapIBD/",cm,"cm/all_",cm,".ibd.global.cut0.088.RData")))
if(cm2 == "5") ibd <- ibd + 0.002 * Diagonal(n = nrow(ibd)) # positive definite
IBD.chol <- chol(ibd)
go_on <- TRUE
while(go_on) {
cat("Generating phenotype data...\n")
age <- rnorm(N, 50, 5)
sex <- rbinom(N, 1, 0.5)
raneff <- as.numeric(crossprod(IBD.chol, rnorm(N)))
ranerr <- rnorm(N)
y <- 0.05*age+0.5*sex+raneff+ranerr
pheno <- data.frame(id = ids, y = y, age = age, sex = sex)
pheno <- pheno[sample(1:N, 400000),]
obj <- glmmkin(y ~ age + sex, data = pheno, id = "id", kins = ibd, family = gaussian(link = "identity"), verbose = F)
if(all(obj$theta>0)) go_on <- FALSE
}
print(obj$theta)
match.idx1 <- match(pheno$id, rownames(ibd))
match.idx2 <- match(pheno$id, colnames(ibd))
print(table(match.idx1 == match.idx2))
tmpibd <- ibd[match.idx1, match.idx2]
IBD.chol <- chol(tmpibd)
Sigma <- obj$theta[1]*Diagonal(nrow(tmpibd))+obj$theta[2]*tmpibd
ncovar <- ncol(obj$X)
obj1 <- glmmkin2randomvec(obj, Z = list(t(IBD.chol)), N.randomvec = N.randomvec)
offset <- crossprod(obj1$random.vectors, crossprod(Sigma, obj1$random.vectors))
rm(obj, tmpibd, IBD.chol, Sigma)
gc()
ibd.file.list <- list.files(paste0("hapIBD/",cm,"cm/chunk.1cM"), pattern = ".RData", full.names = TRUE)
t1 <- proc.time()
out <- ibd_mapping_batch(obj1, offset = offset, ibd.ids = ids, ibd.file.list = ibd.file.list, ncores = ncores)
t2 <- proc.time()
write.table(out, paste0("hapIBD.seed", seed, "_",cm,"cm_run",runs,".txt"), quote=F, row.names=F, col.names=T)
write.table(t(c(t2-t1)), paste0("hapIBD.seed", seed, "_",cm,"cm_run",runs,".time"), quote=F, row.names=F, col.names=T)
