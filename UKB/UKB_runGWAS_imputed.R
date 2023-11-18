Sys.setenv(MKL_NUM_THREADS = 1)
library(GMMAT)
library(Matrix)
library(doMC)

dict <- cbind(c(2:4, 35:37), c(48:50, 20015, 21001, 21002))

cmd <- commandArgs(T)
pheno.idx <- as.numeric(cmd[1]) # 48, 49, 50, 20015, 21001, 21002
cat("Analyzing phenotype", pheno.idx, "...\n")
runs <- 1
ncores <- as.numeric(cmd[2]) # 70
cat("Using", ncores, "cores...\n")
cm2 <- 3

global.ibd.dir <- paste0(cm2,"cm/",runs,"/")
geno.dir <- "data/UKB/imputed/"

rank_norm <- function(x) {
    naidx <- is.na(x)
    xx <- x[!naidx]
    n <- length(xx)
    x[!naidx] <- qnorm((rank(xx) - 0.5)/n)
    return(x)
}

pheno.file <- "data/UKB/pheno/ukb.baseline.WhiteBritish.qced.csv"
#eid,gender,sex,age,age2,age_sex,age2_sex,ethnicity,caucasian,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,PC30,PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40,34.0.0,48.0.0,49.0.0,50.0.0,53.0.0,84.0.0,87.0.0,135.0.0,137.0.0,874.0.0,894.0.0,914.0.0,1070.0.0,1080.0.0,1090.0.0,2966.0.0,2976.0.0,3160.0.0,3627.0.0,3659.0.0,3761.0.0,3786.0.0,3894.0.0,3992.0.0,4012.0.0,4022.0.0,4056.0.0,20001.0.0,20002.0.0,20003.0.0,20006.0.0,20007.0.0,20008.0.0,20009.0.0,20015.0.0,21001.0.0,21002.0.0,22000.0.0,22003.0.0,22004.0.0,22005.0.0,22007.0.0,22008.0.0,22011.0.0,22012.0.0,22013.0.0,22014.0.0,22015.0.0,22182.0.0,40000.0.0,40001.0.0,40002.0.0,40005.0.0,40006.0.0,40007.0.0,40008.0.0,40009.0.0,40010.0.0,40011.0.0,40013.0.0,41066.0.0,41067.0.0,41068.0.0,41069.0.0,41070.0.0,41071.0.0,41072.0.0,41073.0.0,41074.0.0,41075.0.0,41076.0.0,41078.0.0,41079.0.0,41080.0.0,41082.0.0,41083.0.0,41084.0.0,41085.0.0,41086.0.0,41087.0.0,41088.0.0,41089.0.0,41090.0.0,41091.0.0,41092.0.0,41093.0.0,41094.0.0,41095.0.0,41096.0.0,41097.0.0,41098.0.0,41099.0.0,41100.0.0,41101.0.0,41102.0.0,41103.0.0,41104.0.0,41106.0.0,41107.0.0,41108.0.0,41109.0.0,41110.0.0,41111.0.0,41112.0.0,41116.0.0,41117.0.0,41118.0.0,41119.0.0,41120.0.0,41121.0.0,41123.0.0,41125.0.0,41132.0.0,41142.0.0,41143.0.0,41144.0.0,41146.0.0,41148.0.0,41200.0.0,41201.0.0,41202.0.0,41203.0.0,41204.0.0,41205.0.0,41207.0.0,41208.0.0,41210.0.0,41211.0.0,41212.0.0,41213.0.0,41229.0.0,41230.0.0,41233.0.0,41235.0.0,41236.0.0,41237.0.0,41238.0.0,41239.0.0,41240.0.0,41241.0.0,41242.0.0,41243.0.0,41245.0.0,41246.0.0,41248.0.0,41249.0.0,41251.0.0,41252.0.0

ids <- read.table("ids_after_filtering")[,1]
ids_exclude <- read.table("withdrawal.csv")[,1]
ids <- ids[!ids %in% ids_exclude]
N <- length(ids)
print(N)

ibd <- get(load(paste0(global.ibd.dir, "all_",cm2,".max.global.cut0.088.RData")))

i <- dict[dict[,2]==pheno.idx, 1]
pheno <- read.csv(pipe(paste0("cut -d, -f1,3-7,10-49,", 49+i, " ", pheno.file)), as.is = T)
pheno <- subset(pheno, eid %in% ids)
pheno <- pheno[!apply(is.na(pheno), 1, any), ]
pheno.name <- names(pheno)[length(names(pheno))]
print(pheno.name)
mdl <- try(lm(as.formula(paste0(pheno.name, " ~ sex + age + age2 + age_sex + age2_sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),data=pheno))
if(class(mdl) == "try-error") stop("Error: cannot fit the linear model")
pheno$res <- rank_norm(residuals(mdl))
t1 <- proc.time()
obj <- try(glmmkin(res ~ sex + age + age2 + age_sex + age2_sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno, id = "eid", kins = ibd, family = gaussian(link = "identity"), verbose = T))
if(class(obj) == "try-error") stop("Error: cannot fit the LMM")
cat("Number of individuals in glmmkin:", length(obj$id_include), "\n")
cat("Null model: theta =", obj$theta, "\n")
if(obj$theta[1]==0) stop("Error: residual variance component estimate is 0")
for(chr in 1:22) {
print(chr)
system("date")
glmm.score(obj, infile = paste0(geno.dir, "chr",chr,"_v3.impute2.gds"), outfile = paste0("GWAS_imputed.chr",chr,".", pheno.name, ".txt"), MAF.range = c(100/2/length(obj$residuals), 0.5), miss.cutoff = 0.05, ncores = ncores)
}
t2 <- proc.time()
write.table(t(c(t2-t1)), paste0("GWAS_imputed.", pheno.name, ".time"), quote=F, row.names=F, col.names=T)
