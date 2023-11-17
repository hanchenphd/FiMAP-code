# S1 Fig. Quantile-quantile plot of finite-sample and asymptotic IBD mapping variance component test p values under the null hypothesis.
# S2 Fig. Scatter plot of finite-sample and asymptotic IBD mapping variance component test p values under the null hypothesis.

cols <- c("#e69d00", "#56b3e9", "#009e74", "#f0e442", "#0071b2", "#d55c00", "#cc79a7")
# https://doi.org/10.1038/nmeth.1618
# Nat Methods. 2011 Jul;8(7):525.
# Points of view: Avoiding color.
# Wong B, PMID: 21850730

thresh <- 1e-50
seed1 <- 12345 # exact
run1 <- 1
seed2 <- 12345 # B = 100 1000 10000
run2 <- 1

library(tidyverse)
library(ggplot2)
library(cowplot)

outfile1 <- "S1_Fig.eps"
outfile2 <- "S2_Fig.eps"

match.loose <- function(a, b) sapply(a, function(x) which.min(abs(x-b)))
phenovec <- c(48, 49, 50, 20015, 21001, 21002)
ii <- 1
pheno <- phenovec[ii]

label <- "10k"
rapid <- read.table(paste0("3cm/",run1,"/seed",seed1,"_random10000_simIBD_",label,"_exact_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid$p.value[rapid$p.value < thresh] <- thresh
B <- 100
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap$group <- "B = 100"
datap.a <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a$group <- "B = 100"
B <- 1000
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap2 <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap2$group <- "B = 1000"
datap.a2 <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a2$group <- "B = 1000"
datap <- rbind(datap, datap2)
datap.a <- rbind(datap.a, datap.a2)
B <- 10000
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap2 <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap2$group <- "B = 10000"
datap.a2 <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a2$group <- "B = 10000"
datap <- rbind(datap, datap2)
datap.a <- rbind(datap.a, datap.a2)
datap$group <- as.factor(datap$group)
datap.a$group <- as.factor(datap.a$group)

p1_10k <- ggplot(datap, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap$p1, na.rm=T))) +
    ylim(0, -log10(min(datap$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = bquote("FiMAP" ~ ~-log[10](italic(p))), x = bquote("Dense" ~ hat(P) ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6)]) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")

p3_10k <- ggplot(rapid, aes(x = -log10(p.value), y = -log10(p.value.asymptotic))) +
    geom_point() +
    xlim(0, -log10(min(rapid$p.value, na.rm=T))) +
    ylim(0, -log10(min(rapid$p.value.asymptotic, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = bquote("Finite-sample" ~ ~-log[10](italic(p))), y = bquote("Asymptotic" ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[5]) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")

p4_10k <- ggplot(datap.a, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap.a$p1, na.rm=T))) +
    ylim(0, -log10(min(datap.a$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = bquote("FiMAP" ~ ~-log[10](italic(p))), x = bquote("Dense" ~ hat(P) ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6)],name="random matrix") +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")

datap_10k <- datap
p2_10k <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1), las=1)
n0 <- nrow(datap_10k)/3
select0 <- 1:n0
lobs <- -log10(sort(datap_10k$p2[datap_10k$group == "B = 100"]))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(p))), ylab=expression(Observed ~ ~ -log[10](italic(p))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.33)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(datap_10k$p2[datap_10k$group == "B = 1000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(datap_10k$p2[datap_10k$group == "B = 10000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
lobs <- -log10(sort(datap_10k$p1[datap_10k$group == "B = 100"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[5])
}

datap_10k.a <- datap.a
p5_10k <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1), las=1)
n0 <- nrow(datap_10k.a)/3
select0 <- 1:n0
lobs <- -log10(sort(datap_10k.a$p2[datap_10k.a$group == "B = 100"]))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(p))), ylab=expression(Observed ~ ~ -log[10](italic(p))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.33)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(datap_10k.a$p2[datap_10k.a$group == "B = 1000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(datap_10k.a$p2[datap_10k.a$group == "B = 10000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
lobs <- -log10(sort(datap_10k.a$p1[datap_10k.a$group == "B = 100"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[5])
}

label <- "100k"
rapid <- read.table(paste0("3cm/",run1,"/seed",seed1,"_random10000_simIBD_",label,"_exact_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid$p.value[rapid$p.value < thresh] <- thresh
B <- 100
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap$group <- "B = 100"
datap.a <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a$group <- "B = 100"
B <- 1000
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap2 <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap2$group <- "B = 1000"
datap.a2 <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a2$group <- "B = 1000"
datap <- rbind(datap, datap2)
datap.a <- rbind(datap.a, datap.a2)
B <- 10000
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap2 <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap2$group <- "B = 10000"
datap.a2 <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a2$group <- "B = 10000"
datap <- rbind(datap, datap2)
datap.a <- rbind(datap.a, datap.a2)
datap$group <- as.factor(datap$group)
datap.a$group <- as.factor(datap.a$group)

p1_100k <- ggplot(datap, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap$p1, na.rm=T))) +
    ylim(0, -log10(min(datap$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = bquote("FiMAP" ~ ~-log[10](italic(p))), x = bquote("Dense" ~ hat(P) ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6)]) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")

p3_100k <- ggplot(rapid, aes(x = -log10(p.value), y = -log10(p.value.asymptotic))) +
    geom_point() +
    xlim(0, -log10(min(rapid$p.value, na.rm=T))) +
    ylim(0, -log10(min(rapid$p.value.asymptotic, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = bquote("Finite-sample" ~ ~-log[10](italic(p))), y = bquote("Asymptotic" ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[5]) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")

p4_100k <- ggplot(datap.a, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap.a$p1, na.rm=T))) +
    ylim(0, -log10(min(datap.a$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = bquote("FiMAP" ~ ~-log[10](italic(p))), x = bquote("Dense" ~ hat(P) ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6)],name="random matrix") +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")

datap_100k <- datap
p2_100k <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1), las=1)
n0 <- nrow(datap_100k)/3
select0 <- 1:n0
lobs <- -log10(sort(datap_100k$p2[datap_100k$group == "B = 100"]))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(p))), ylab=expression(Observed ~ ~ -log[10](italic(p))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.33)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(datap_100k$p2[datap_100k$group == "B = 1000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(datap_100k$p2[datap_100k$group == "B = 10000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
lobs <- -log10(sort(datap_100k$p1[datap_100k$group == "B = 100"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[5])
}

datap_100k.a <- datap.a
p5_100k <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1), las=1)
n0 <- nrow(datap_100k.a)/3
select0 <- 1:n0
lobs <- -log10(sort(datap_100k.a$p2[datap_100k.a$group == "B = 100"]))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(p))), ylab=expression(Observed ~ ~ -log[10](italic(p))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.33)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(datap_100k.a$p2[datap_100k.a$group == "B = 1000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(datap_100k.a$p2[datap_100k.a$group == "B = 10000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
lobs <- -log10(sort(datap_100k.a$p1[datap_100k.a$group == "B = 100"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[5])
}

label <- "1m"
rapid <- read.table(paste0("3cm/",run1,"/seed",seed1,"_random10000_simIBD_",label,"_exact_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid$p.value[rapid$p.value < thresh] <- thresh
B <- 100
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap$group <- "B = 100"
datap.a <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a$group <- "B = 100"
B <- 1000
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap2 <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap2$group <- "B = 1000"
datap.a2 <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a2$group <- "B = 1000"
datap <- rbind(datap, datap2)
datap.a <- rbind(datap.a, datap.a2)
B <- 10000
rapid2 <- read.table(paste0("3cm/",run2,"/seed",seed2,"_random10000_simIBD_",label,"_B",B,"_X", pheno, ".0.0.txt"), header = T, as.is = T)
rapid2$p.value[rapid2$p.value < thresh] <- thresh
datap2 <- data.frame(p1 = rapid$p.value, p2 = rapid2$p.value)
datap2$group <- "B = 10000"
datap.a2 <- data.frame(p1 = rapid$p.value.asymptotic, p2 = rapid2$p.value.asymptotic)
datap.a2$group <- "B = 10000"
datap <- rbind(datap, datap2)
datap.a <- rbind(datap.a, datap.a2)
datap$group <- as.factor(datap$group)
datap.a$group <- as.factor(datap.a$group)

p1 <- ggplot(datap, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap$p1, na.rm=T))) +
    ylim(0, -log10(min(datap$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = bquote("FiMAP" ~ ~-log[10](italic(p))), x = bquote("Dense" ~ hat(P) ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6)],name="random matrix") +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size=4)))
legend1 <- get_legend(p1)
p1_1m <- p1 + theme(legend.position = "none")

p3_1m <- ggplot(rapid, aes(x = -log10(p.value), y = -log10(p.value.asymptotic))) +
    geom_point() +
    xlim(0, -log10(min(rapid$p.value, na.rm=T))) +
    ylim(0, -log10(min(rapid$p.value.asymptotic, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = bquote("Finite-sample" ~ ~-log[10](italic(p))), y = bquote("Asymptotic" ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[5]) +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.position = "none")

p4 <- ggplot(datap.a, aes(x = -log10(p1), y = -log10(p2), color = group)) +
    geom_point() +
    xlim(0, -log10(min(datap.a$p1, na.rm=T))) +
    ylim(0, -log10(min(datap.a$p2, na.rm=T))) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = bquote("FiMAP" ~ ~-log[10](italic(p))), x = bquote("Dense" ~ hat(P) ~ ~-log[10](italic(p)))) +
    scale_color_manual(values=cols[c(1,2,6)],name="random matrix") +
    theme(axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    plot.margin = margin(50, 15, 0, 20),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size=4)))
legend4 <- get_legend(p4)
p4_1m <- p4 + theme(legend.position = "none")

p2_1m <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1), las=1)
n0 <- nrow(datap)/3
select0 <- 1:n0
lobs <- -log10(sort(datap$p2[datap$group == "B = 100"]))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(p))), ylab=expression(Observed ~ ~ -log[10](italic(p))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.33)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(datap$p2[datap$group == "B = 1000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(datap$p2[datap$group == "B = 10000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
lobs <- -log10(sort(datap$p1[datap$group == "B = 100"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[5])
legend("topleft", legend=c(bquote(B == 100), bquote(B == 1000), bquote(B == 10000), bquote("Dense" ~ hat(P))), lty=1, lwd=3, col=cols[c(1,2,6,5)], bty="n", cex=1.2)
}

p5_1m <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1), las=1)
n0 <- nrow(datap.a)/3
select0 <- 1:n0
lobs <- -log10(sort(datap.a$p2[datap.a$group == "B = 100"]))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(p))), ylab=expression(Observed ~ ~ -log[10](italic(p))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.33)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(datap.a$p2[datap.a$group == "B = 1000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(datap.a$p2[datap.a$group == "B = 10000"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
lobs <- -log10(sort(datap.a$p1[datap.a$group == "B = 100"]))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[5])
legend("topleft", legend=c(bquote(B == 100), bquote(B == 1000), bquote(B == 10000), bquote("Dense" ~ hat(P))), lty=1, lwd=3, col=cols[c(1,2,6,5)], bty="n", cex=1.2)
}


row1 <- plot_grid(p1_10k, p1_100k, p1_1m, legend1, labels = c("A","B","C"), rel_widths = c(1,1,1,0.4), nrow=1, label_size = 30)
row2 <- plot_grid(p2_10k, p2_100k, p2_1m, labels = c("A","B","C"), nrow=1, label_size = 30)
row3 <- plot_grid(p3_10k, p3_100k, p3_1m, labels = c("D","E","F"), nrow=1, label_size = 30)
row4 <- plot_grid(p4_10k, p4_100k, p4_1m, legend4, labels = c("G","H","I"), rel_widths = c(1,1,1,0.4), nrow=1, label_size = 30)
row5 <- plot_grid(p5_10k, p5_100k, p5_1m, labels = c("D","E","F"), nrow=1, label_size = 30)
plot_grid(row2, row5, ncol = 1)
ggsave(outfile1, width = 18, height = 12, bg = "white")
plot_grid(row1, row3, row4, ncol = 1)
ggsave(outfile2, width = 18, height = 18, bg = "white")

