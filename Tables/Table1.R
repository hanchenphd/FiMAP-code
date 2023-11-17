# Table 1. Proportions of power simulation replicates in which both FiMAP and GWAS or only one approach identified the association.
rsq <- 0.02
alpha <- 0.05/3403
alpha.gwas <- 5e-8

data <- NULL
for(seed in 1:40) {
infile1 <- paste0("untyped/rsq_",rsq,"_seed",seed,".txt")
tmp <- read.table(infile1, header=T)
data <- rbind(data, tmp)
}
data1 <- data[!is.na(data$p.svt.min1),c("p.value3_1", "p.value.hapIBD3_1", "p.value5_1", "p.value.hapIBD5_1", "p.value10_1", "p.value.hapIBD10_1", "p.svt.min1", "p.svt.cauchy1", "p.svt.imputed.min1", "p.svt.imputed.cauchy1")]
colnames(data1) <- c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10", "p.svt.min", "p.svt.cauchy", "p.svt.imputed.min", "p.svt.imputed.cauchy")
data2 <- data[!is.na(data$p.svt.min2),c("p.value3_2", "p.value.hapIBD3_2", "p.value5_2", "p.value.hapIBD5_2", "p.value10_2", "p.value.hapIBD10_2", "p.svt.min2", "p.svt.cauchy2", "p.svt.imputed.min2", "p.svt.imputed.cauchy2")]
colnames(data2) <- c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10", "p.svt.min", "p.svt.cauchy", "p.svt.imputed.min", "p.svt.imputed.cauchy")
res <- data.frame(Test = c("Both", "FiMAP only", "GWAS only"), Power = c(sum(data$p.value3<alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data), sum(data$p.value3<alpha & data$p.svt.imputed.min>=alpha.gwas)/nrow(data), sum(data$p.value3>=alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data)))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))
print(res) # causal

data <-rbind(data1, data2)
res <- data.frame(Test = c("Both", "FiMAP only", "GWAS only"), Power = c(sum(data$p.value3<alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data), sum(data$p.value3<alpha & data$p.svt.imputed.min>=alpha.gwas)/nrow(data), sum(data$p.value3>=alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data)))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))
print(res) # neighboring

data <- NULL
for(seed in 1:40) {
infile2 <- paste0("haplotype/rsq_",rsq,"_seed",seed,".txt")
tmp <- read.table(infile2, header=T)
data <- rbind(data, tmp)
}
data1 <- data[!is.na(data$p.svt.min1),c("p.value3_1", "p.value.hapIBD3_1", "p.value5_1", "p.value.hapIBD5_1", "p.value10_1", "p.value.hapIBD10_1", "p.svt.min1", "p.svt.cauchy1", "p.svt.imputed.min1", "p.svt.imputed.cauchy1")]
colnames(data1) <- c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10", "p.svt.min", "p.svt.cauchy", "p.svt.imputed.min", "p.svt.imputed.cauchy")
data2 <- data[!is.na(data$p.svt.min2),c("p.value3_2", "p.value.hapIBD3_2", "p.value5_2", "p.value.hapIBD5_2", "p.value10_2", "p.value.hapIBD10_2", "p.svt.min2", "p.svt.cauchy2", "p.svt.imputed.min2", "p.svt.imputed.cauchy2")]
colnames(data2) <- c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10", "p.svt.min", "p.svt.cauchy", "p.svt.imputed.min", "p.svt.imputed.cauchy")
res <- data.frame(Test = c("Both", "FiMAP only", "GWAS only"), Power = c(sum(data$p.value3<alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data), sum(data$p.value3<alpha & data$p.svt.imputed.min>=alpha.gwas)/nrow(data), sum(data$p.value3>=alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data)))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))
print(res) # causal

data <-rbind(data1, data2)
res <- data.frame(Test = c("Both", "FiMAP only", "GWAS only"), Power = c(sum(data$p.value3<alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data), sum(data$p.value3<alpha & data$p.svt.imputed.min>=alpha.gwas)/nrow(data), sum(data$p.value3>=alpha & data$p.svt.imputed.min<alpha.gwas)/nrow(data)))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))
print(res) # neighboring
