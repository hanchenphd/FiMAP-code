# Fig 3. Quantile-quantile plot of finite-sample FiMAP p values under the null hypothesis.

library(ggplot2)
library(cowplot)
match.loose <- function(a, b) sapply(a, function(x) which.min(abs(x-b)))
cols <- c("#e69d00", "#56b3e9", "#009e74", "#f0e442", "#0071b2", "#d55c00", "#cc79a7")
# https://doi.org/10.1038/nmeth.1618
# Nat Methods. 2011 Jul;8(7):525.
# Points of view: Avoiding color.
# Wong B, PMID: 21850730

outfile <- "Fig3.eps"

infiles3 <- paste0("seed",1:50,"_3cm_run1.txt")
infiles5 <- paste0("seed",1:50,"_5cm_run1.txt")
infiles10 <- paste0("seed",1:50,"_10cm_run1.txt")

data3 <- NULL
for(f in infiles3) data3 <- rbind(data3, read.table(f, header=T))
data5 <- NULL
for(f in infiles5) data5 <- rbind(data5, read.table(f, header=T))
data10 <- NULL
for(f in infiles10) data10 <- rbind(data10, read.table(f, header=T))
print(nrow(data3));print(nrow(data5));print(nrow(data10))
lambda3.as <- qchisq(median(data3$p.value.asymptotic), 1, lower = F)/qchisq(0.5, 1)
lambda3.fs <- qchisq(median(data3$p.value), 1, lower = F)/qchisq(0.5, 1)
lambda5.as <- qchisq(median(data5$p.value.asymptotic), 1, lower = F)/qchisq(0.5, 1)
lambda5.fs <- qchisq(median(data5$p.value), 1, lower = F)/qchisq(0.5, 1)
lambda10.as <- qchisq(median(data10$p.value.asymptotic), 1, lower = F)/qchisq(0.5, 1)
lambda10.fs <- qchisq(median(data10$p.value), 1, lower = F)/qchisq(0.5, 1)
print(lambda3.as); print(lambda3.fs); print(lambda5.as); print(lambda5.fs); print(lambda10.as); print(lambda10.fs)

p1 <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1))
n0 <- nrow(data3)
select0 <- 1:n0
lobs <- -log10(sort(data3$p.value))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(P))), ylab=expression(Observed ~ ~ -log[10](italic(P))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=2, cex.axis=1.8)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(data5$p.value))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(data10$p.value))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
legend("topleft", legend=c(as.expression(substitute(paste("3cM, ", lambda[GC]==ll, sep=""), list(ll=round(lambda3.fs,3)))),as.expression(substitute(paste("5cM, ", lambda[GC]==ll, sep=""), list(ll=round(lambda5.fs,3)))),as.expression(substitute(paste("10cM, ", lambda[GC]==ll, sep=""), list(ll=round(lambda10.fs,3))))), lty=1, lwd=3, col=cols[c(1,2,6)], bty="n", cex=1.2)
}

infiles3 <- paste0("hapIBD.seed",1:50,"_3cm_run1.txt")
infiles5 <- paste0("hapIBD.seed",1:50,"_5cm_run1.txt")
infiles10 <- paste0("hapIBD.seed",1:50,"_10cm_run1.txt")

hapIBD.data3 <- NULL
for(f in infiles3) hapIBD.data3 <- rbind(hapIBD.data3, read.table(f, header=T))
hapIBD.data5 <- NULL
for(f in infiles5) hapIBD.data5 <- rbind(hapIBD.data5, read.table(f, header=T))
hapIBD.data10 <- NULL
for(f in infiles10) hapIBD.data10 <- rbind(hapIBD.data10, read.table(f, header=T))
print(nrow(hapIBD.data3));print(nrow(hapIBD.data5));print(nrow(hapIBD.data10))
hapIBD.lambda3.as <- qchisq(median(hapIBD.data3$p.value.asymptotic), 1, lower = F)/qchisq(0.5, 1)
hapIBD.lambda3.fs <- qchisq(median(hapIBD.data3$p.value), 1, lower = F)/qchisq(0.5, 1)
hapIBD.lambda5.as <- qchisq(median(hapIBD.data5$p.value.asymptotic), 1, lower = F)/qchisq(0.5, 1)
hapIBD.lambda5.fs <- qchisq(median(hapIBD.data5$p.value), 1, lower = F)/qchisq(0.5, 1)
hapIBD.lambda10.as <- qchisq(median(hapIBD.data10$p.value.asymptotic), 1, lower = F)/qchisq(0.5, 1)
hapIBD.lambda10.fs <- qchisq(median(hapIBD.data10$p.value), 1, lower = F)/qchisq(0.5, 1)
print(hapIBD.lambda3.as); print(hapIBD.lambda3.fs); print(hapIBD.lambda5.as); print(hapIBD.lambda5.fs); print(hapIBD.lambda10.as); print(hapIBD.lambda10.fs)

p2 <- function() {
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.5,0.1))
n0 <- nrow(hapIBD.data3)
select0 <- 1:n0
lobs <- -log10(sort(hapIBD.data3$p.value))
funnelx <- seq(log10(1+n0), 0, by=-0.01)
Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](italic(P))), ylab=expression(Observed ~ ~ -log[10](italic(P))), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=2, cex.axis=1.8)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.5, col=cols[1])
lobs <- -log10(sort(hapIBD.data5$p.value))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[2])
lobs <- -log10(sort(hapIBD.data10$p.value))
lobs <- lobs[Select0]
lexp <- -log10(select0[Select0]/(1+n0))
points(lexp, lobs, pch=16, cex=0.5, col=cols[6])
legend("topleft", legend=c(as.expression(substitute(paste("3cM, ", lambda[GC]==ll, sep=""), list(ll=round(hapIBD.lambda3.fs,3)))),as.expression(substitute(paste("5cM, ", lambda[GC]==ll, sep=""), list(ll=round(hapIBD.lambda5.fs,3)))),as.expression(substitute(paste("10cM, ", lambda[GC]==ll, sep=""), list(ll=round(hapIBD.lambda10.fs,3))))), lty=1, lwd=3, col=cols[c(1,2,6)], bty="n", cex=1.2)
}

plot_grid(p1, p2, ncol = 2, align = "hv", labels = "AUTO", label_size = 30)
ggsave(outfile, width=12, height=6, bg="white")
