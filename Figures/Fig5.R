# Fig 5. Comparison of finite-sample FiMAP p values from RaPID IBD segments before and after conditioning on GWAS tag variants in the testing window and flanking regions.

cols <- c("#e69d00", "#56b3e9", "#009e74", "#f0e442", "#0071b2", "#d55c00", "#cc79a7")
# https://doi.org/10.1038/nmeth.1618
# Nat Methods. 2011 Jul;8(7):525.
# Points of view: Avoiding color.
# Wong B, PMID: 21850730

thresh <- 1e-50
seed1 <- 12345

library(tidyverse)
library(ggplot2)
library(cowplot)

alpha <- "1e-06"
outfile <- "Fig5.eps"

ps <- list()
phenovec <- c(48, 49, 50, 20015, 21001, 21002)
for(ii in 1:6) {
pheno <- phenovec[ii]
rapid <- read.table(paste0("3cm/1/seed",seed1,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid$p.value[rapid$p.value < thresh] <- thresh
rapid2 <- read.table(paste0("3cm/1/conditional_imputed_GDS_pheno",pheno,"_maf1e-04_missrate0.05_qual0.3_flank3cM_alpha",alpha,".out"), header = T, as.is = T,sep="\t")
rapid2$p.value.unconditional[rapid2$p.value.unconditional < thresh] <- thresh
rapid2$p.value.conditional[rapid2$p.value.conditional < thresh] <- thresh
rapid2$group <- ifelse(rapid2$p.value.unconditional<0.05/nrow(rapid), ifelse(rapid2$p.value.conditional<0.05/nrow(rapid), 6, 5), NA)

tmp <- ggplot(rapid2, aes(x = -log10(p.value.unconditional), y = -log10(p.value.conditional), color = as.factor(group))) +
    geom_point() +
    xlim(0, -log10(min(rapid2$p.value.unconditional, na.rm=T))) +
    ylim(0, -log10(min(rapid2$p.value.conditional, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = -log10(0.05/nrow(rapid)), linetype = 2) +
    geom_hline(yintercept = -log10(0.05/nrow(rapid)), linetype = 2) +
    labs(x = expression(Unconditional ~ ~-log[10](italic(p))), y = expression(Conditional ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[5:6]) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(40, 10, 0, 20),
    legend.position = "none")
ps <- c(ps, list(tmp))
}

plot_grid(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ps[[6]], ncol = 3, align = "hv", labels = "AUTO", label_size = 30)
ggsave(outfile, width=12, height=7.5)
