# Fig 4. Power comparison of FiMAP and GWAS.

rsq <- 0.02
outfile <- "Fig4.eps"
alpha <- 0.05/3403
alpha.gwas <- 5e-8

library(tidyverse)
library(ggplot2)
library(cowplot)

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
res <- data.frame(Test = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM"), Power = colSums(data[,c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10")]<alpha)/nrow(data))
res <- rbind(res, data.frame(Test = c("GWAS"), Power = c(sum(data$p.svt.imputed.min<alpha.gwas)/nrow(data))))
res$Test <- factor(res$Test, levels = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM","GWAS"))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))

# causal
p1 <- ggplot(res, aes(x = Test, y = Power, fill = Test)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Power-se, ymax = Power+se), width = 0.5) +
    ylim(0, 1) +
    theme(axis.title = element_text(size=20),
    axis.text.x = element_text(angle=60,hjust=1,size=16),
    axis.text.y = element_text(size=16),
    plot.margin = margin(20, 5, 0, 20),
    legend.position = "none")

data <-rbind(data1, data2)
res <- data.frame(Test = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM"), Power = colSums(data[,c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10")]<alpha)/nrow(data))
res <- rbind(res, data.frame(Test = c("GWAS"), Power = c(sum(data$p.svt.imputed.min<alpha.gwas)/nrow(data))))
res$Test <- factor(res$Test, levels = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM","GWAS"))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))

# neighboring
p2 <- ggplot(res, aes(x = Test, y = Power, fill = Test)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Power-se, ymax = Power+se), width = 0.5) +
    ylim(0, 1) +
    theme(axis.title = element_text(size=20),
    axis.text.x = element_text(angle=60,hjust=1,size=16),
    axis.text.y = element_text(size=16),
    plot.margin = margin(20, 5, 0, 20),
    legend.position = "none")

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
res <- data.frame(Test = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM"), Power = colSums(data[,c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10")]<alpha)/nrow(data))
res <- rbind(res, data.frame(Test = c("GWAS"), Power = c(sum(data$p.svt.imputed.min<alpha.gwas)/nrow(data))))
res$Test <- factor(res$Test, levels = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM","GWAS"))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))

# causal
p3 <- ggplot(res, aes(x = Test, y = Power, fill = Test)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Power-se, ymax = Power+se), width = 0.5) +
    ylim(0, 1) +
    theme(axis.title = element_text(size=20),
    axis.text.x = element_text(angle=60,hjust=1,size=16),
    axis.text.y = element_text(size=16),
    plot.margin = margin(20, 5, 0, 20),
    legend.position = "none")

data <-rbind(data1, data2)
res <- data.frame(Test = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM"), Power = colSums(data[,c("p.value3", "p.value.hapIBD3", "p.value5", "p.value.hapIBD5", "p.value10", "p.value.hapIBD10")]<alpha)/nrow(data))
res <- rbind(res, data.frame(Test = c("GWAS"), Power = c(sum(data$p.svt.imputed.min<alpha.gwas)/nrow(data))))
res$Test <- factor(res$Test, levels = c("RaPID 3cM","hap-IBD 3cM","RaPID 5cM","hap-IBD 5cM","RaPID 10cM","hap-IBD 10cM","GWAS"))
res$se <- sqrt(res$Power * (1 - res$Power)/nrow(data))

# neighboring
p4 <- ggplot(res, aes(x = Test, y = Power, fill = Test)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Power-se, ymax = Power+se), width = 0.5) +
    ylim(0, 1) +
    theme(axis.title = element_text(size=20),
    axis.text.x = element_text(angle=60,hjust=1,size=16),
    axis.text.y = element_text(size=16),
    plot.margin = margin(20, 5, 0, 20),
    legend.position = "none")

plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", labels = "AUTO", label_size = 30)
ggsave(outfile, width=12, height=12)
