# S4 Fig. Comparison of finite-sample FiMAP p values using two separate RaPID IBD segment calls with different random number seeds.
# S5 Fig. Comparison of finite-sample FiMAP p values using two random matrices with different random number seeds.

cols <- c("#e69d00", "#56b3e9", "#009e74", "#f0e442", "#0071b2", "#d55c00", "#cc79a7")
# https://doi.org/10.1038/nmeth.1618
# Nat Methods. 2011 Jul;8(7):525.
# Points of view: Avoiding color.
# Wong B, PMID: 21850730

thresh <- 1e-50
seed1 <- 12345
run1 <- 1
seed2 <- 12345
run2 <- 2

library(tidyverse)
library(ggplot2)
library(cowplot)

outfile <- "S4_Fig.eps"

ps <- list()
phenovec <- c(48, 49, 50, 20015, 21001, 21002)
for(ii in 1:6) {
pheno <- phenovec[ii]
rapid <- read.table(paste0("3cm/",run1,"/seed",seed1,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid$p.value[rapid$p.value < thresh] <- thresh
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap$group <- factor(ifelse(datap$p1>=0.05/nrow(datap) & datap$p2>=0.05/nrow(datap),1,ifelse(datap$p1>=0.05/nrow(datap),2,ifelse(datap$p2>=0.05/nrow(datap),5,6))),levels=c(1,2,6,5))

tmp <- ggplot(datap, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap$p1, na.rm=T))) +
    ylim(0, -log10(min(datap$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = -log10(0.05/nrow(datap)), linetype = 2) +
    geom_hline(yintercept = -log10(0.05/nrow(datap)), linetype = 2) +
    labs(x = expression(RaPID~ ~ seed~ ~ 1 ~ ~-log[10](italic(p))), y = expression(RaPID~ ~ seed~ ~ 2 ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6,5)]) +
    geom_text(aes(x=0,y=0.98*(-log10(min(p2, na.rm=T))),label=paste("rho","==", round(cor(p1, p2, method = "spearman"),3))),hjust=0,vjust=1,size=6,parse=T,inherit.aes=F) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")
ps <- c(ps, list(tmp))
}

plot_grid(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ps[[6]], ncol = 3, align = "hv", labels = "AUTO", label_size = 30)
ggsave(outfile, width=12, height=7.5)

seed2 <- 54321
run2 <- 1
outfile <- "S5_Fig.eps"

ps <- list()
for(ii in 1:6) {
pheno <- phenovec[ii]
rapid <- read.table(paste0("3cm/",run1,"/seed",seed1,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid$p.value[rapid$p.value < thresh] <- thresh
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap$group <- factor(ifelse(datap$p1>=0.05/nrow(datap) & datap$p2>=0.05/nrow(datap),1,ifelse(datap$p1>=0.05/nrow(datap),2,ifelse(datap$p2>=0.05/nrow(datap),5,6))),levels=c(1,2,6,5))

tmp <- ggplot(datap, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap$p1, na.rm=T))) +
    ylim(0, -log10(min(datap$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = -log10(0.05/nrow(datap)), linetype = 2) +
    geom_hline(yintercept = -log10(0.05/nrow(datap)), linetype = 2) +
    labs(x = expression(FiMAP~ ~ seed~ ~ 1 ~ ~-log[10](italic(p))), y = expression(FiMAP~ ~ seed~ ~ 2 ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6,5)]) +
    geom_text(aes(x=0,y=0.98*(-log10(min(p2, na.rm=T))),label=paste("rho","==", round(cor(p1, p2, method = "spearman"),3))),hjust=0,vjust=1,size=6,parse=T,inherit.aes=F) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")
ps <- c(ps, list(tmp))
}

plot_grid(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ps[[6]], ncol = 3, align = "hv", labels = "AUTO", label_size = 30)
ggsave(outfile, width=12, height=7.5)
