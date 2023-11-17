# S3 Fig. GWAS and finite-sample FiMAP p values from RaPID IBD segments on chromosome 16.

cols <- c("#e69d00", "#56b3e9", "#009e74", "#f0e442", "#0071b2", "#d55c00", "#cc79a7")
# https://doi.org/10.1038/nmeth.1618
# Nat Methods. 2011 Jul;8(7):525.
# Points of view: Avoiding color.
# Wong B, PMID: 21850730

CHR <- "16"
pval.cutoff1 <- 1e-25 # other phenotypes
pval.cutoff2 <- 1e-50 # Standing height

## Code credit: Matthew Flickinger
## https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
library(lattice)
library(ggplot2)
library(cowplot)
manhattan.plot<-function(chr, pos, pvalue, sig.level=NA, annotate=NULL, ann.default=list(), should.thin=T, thin.pos.places=2, thin.logp.places=2, xlab="Chromosome", ylab=expression(-log[10](p-value)), col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
      if (length(chr)==0) stop("chromosome vector is empty")
      if (length(pos)==0) stop("position vector is empty")
      if (length(pvalue)==0) stop("pvalue vector is empty")
      #make sure we have an ordered factor
      if(!is.ordered(chr)) {
	    chr <- ordered(chr)
      } else {
	    chr <- chr[,drop=T]
      }
      #make sure positions are in kbp
      if (any(pos>1e6)) pos<-pos/1e6;
      #calculate absolute genomic position
      #from relative chromosomal positions
      posmin <- tapply(pos,chr, min);
      posmax <- tapply(pos,chr, max);
      posshift <- head(c(0,cumsum(posmax)),-1);
      names(posshift) <- levels(chr)
      genpos <- pos + posshift[chr];
      getGenPos<-function(cchr, cpos) {
     	    p<-posshift[as.character(cchr)]+cpos
	    return(p)
      }
      #parse annotations
      grp <- NULL
      ann.settings <- list()
      label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, col=NULL, fontface=NULL, fontsize=NULL, show=F)
      parse.label<-function(rawval, groupname) {
 	    r<-list(text=groupname)
	    if(is.logical(rawval)) {
		  if(!rawval) {r$show <- F}
	    } else if(is.character(rawval) || is.expression(rawval)) {
		  if(nchar(rawval)>=1) {
			r$text <- rawval
		  }
	    } else if(is.list(rawval)) {
		  r <- modifyList(r, rawval)
	    }
	    return(r)
      }

      if(!is.null(annotate)) {
	    if(is.list(annotate)) {
		  grp <- annotate[[1]]
	    } else {
		  grp <- annotate
	    } 
	    if(!is.factor(grp)) {
		  grp <- factor(grp)
	    }
      }	else {
	    grp <- factor(rep(1, times=length(pvalue)))
      }
  
      ann.settings<-vector("list", length(levels(grp)))
      ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

      if (length(ann.settings)>1) { 
            lcols<-trellis.par.get("superpose.symbol")$col 
     	    lfills<-trellis.par.get("superpose.symbol")$fill
	    for(i in 2:length(levels(grp))) {
		  ann.settings[[i]]<-list(pch=pch, col=lcols[(i-2) %% length(lcols) +1 ], fill=lfills[(i-2) %% length(lfills) +1 ], cex=cex, label=label.default);
		  ann.settings[[i]]$label$show <- T
	    }
	    names(ann.settings)<-levels(grp)
      }
      for(i in 1:length(ann.settings)) {
     	    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
     	    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, parse.label(ann.settings[[i]]$label, levels(grp)[i]))
      }
      if(is.list(annotate) && length(annotate)>1) {
     	    user.cols <- 2:length(annotate)
	    ann.cols <- c()
	    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
     		  ann.cols<-match(names(annotate)[-1], names(ann.settings))
	    } else {
	          ann.cols<-user.cols-1
	    }
	    for(i in seq_along(user.cols)) {
		  if(!is.null(annotate[[user.cols[i]]]$label)) {
			annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, levels(grp)[ann.cols[i]])
		  }
		  ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], annotate[[user.cols[i]]])
	    }
      }
      rm(annotate)

      #reduce number of points plotted
      if(should.thin) {
	    thinned <- unique(data.frame(logp=round(-log10(pvalue),thin.logp.places), pos=round(genpos,thin.pos.places), chr=chr, grp=grp))
	    logp <- thinned$logp
	    genpos <- thinned$pos
	    chr <- thinned$chr
	    grp <- thinned$grp
	    rm(thinned)
      } else {
	    logp <- -log10(pvalue)
      }
      rm(pos, pvalue)
      gc()

      #custom axis to print chromosome names
      axis.chr <- function(side,...) {
	    if(side=="bottom") {
#		  panel.axis(side=side, outside=T, at=((posmax+posmin)/2+posshift), labels=levels(chr), ticks=F, rot=0, check.overlap=F)
		  panel.axis(side=side, outside=T, at=((posmax+posmin)/2+posshift), labels=paste("Chromosome",CHR), text.cex = 1.8, ticks=F, rot=0, check.overlap=F)
	    } else if(side=="top" || side=="right") {
		  panel.axis(side=side, draw.labels=F, ticks=F);
	    } else {
		  axis.default(side=side, ...);
	    }
      }

      #make sure the y-lim covers the range (plus a bit more to look nice)
      prepanel.chr<-function(x,y,...) { 
	    A<-list();
	    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
	    A$ylim=c(0,maxy);
	    A;
      }

#      xyplot(logp~genpos, chr=chr, groups=grp, axis=axis.chr, ann.settings=ann.settings, prepanel=prepanel.chr, scales=list(axs="i"), panel=function(x, y, ..., getgenpos) {
      xyplot(logp~genpos, chr=chr, groups=grp, axis=axis.chr, ann.settings=ann.settings, prepanel=prepanel.chr, scales=list(axs="i",y=list(cex=1.8)), panel=function(x, y, ..., getgenpos) {
	    if(!is.na(sig.level)) {
		  #add significance line (if requested)
		  panel.abline(h=-log10(sig.level), lty=2);
	    }
	    panel.superpose(x, y, ..., getgenpos=getgenpos);
	    if(!is.null(panel.extra)) {
		  panel.extra(x,y, getgenpos, ...)
	    }
      }, panel.groups = function(x,y,..., subscripts, group.number) {
	    A<-list(...)
	    #allow for different annotation settings
	    gs <- ann.settings[[group.number]]
	    A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
	    A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
	    A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
	    A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
	    A$x <- x
	    A$y <- y
	    do.call("panel.xyplot", A)
	    #draw labels (if requested)
	    if(gs$label$show) {
		  gt<-gs$label
		  names(gt)[which(names(gt)=="text")]<-"labels"
		  gt$show<-NULL
		  if(is.character(gt$x) | is.character(gt$y)) {
			peak = which.max(y)
			center = mean(range(x))
			if(is.character(gt$x)) {
			      if(gt$x=="peak") {gt$x<-x[peak]}
			      if(gt$x=="center") {gt$x<-center}
			}
			if(is.character(gt$y)) {
			      if(gt$y=="peak") {gt$y<-y[peak]}
			}
		  }
		  if(is.list(gt$x)) {
			gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
		  }
		  do.call("panel.text", gt)
	    }
      }, xlab=xlab, ylab=ylab, panel.extra=panel.extra, getgenpos=getGenPos, ...);
}

##########################
# Waist circumference (48.0.0)
# Hip circumference (49.0.0)
# Standing height (50.0.0)
# Sitting height (20015.0.0)
# Body mass index (21001.0.0)
# Weight (21002.0.0)

alpha <- "1e-06"
maf <- 0.0001
missrate <- 0.05
qual <- 0.3

infile <- paste0("GWAS_imputed.chr",CHR,".X48.0.0.QUALgt0.3.txt")
data <- read.table(infile, header=T, as.is=T)
data <- subset(data, AF >= maf & AF <= 1 - maf)
data <- subset(data, MISSRATE < missrate)
data <- subset(data, QUAL > qual) # imputed
data <- subset(data, !is.na(PVAL))
data <- data[,c("CHR", "cM", "PVAL")]
colnames(data) <- c("chr", "pos", "p.value")
data$source <- 1

infile <- "3cm/1/seed12345_X48.0.0.txt"
data2 <- read.table(infile, header=T)
data2 <- subset(data2, chr == CHR)
data2$pos <- (data2$start+data2$end)/2
data2 <- data2[, c("chr", "pos", "p.value")]
data2$source <- 2

infile <- paste0("3cm/1/conditional_imputed_GDS_pheno48_maf1e-04_missrate0.05_qual0.3_flank3cM_alpha",alpha,".out")
data3 <- read.table(infile, header=T, as.is=T, sep="\t")
data3 <- subset(data3, chr == CHR)
data3$pos <- (data3$start+data3$end)/2
data3 <- data3[, c("chr", "pos", "p.value.conditional")]
colnames(data3)[3] <- "p.value"
if(nrow(data3)>0) data3$source <- 3

data <- rbind(data, data2, data3)
data$p.value[data$p.value < pval.cutoff1] <- pval.cutoff1
data <- data[order(data$chr, data$pos),]
ann <- factor(data$source, levels = 1:3, labels = c("", "tmp1", "tmp2"))
p1 <- manhattan.plot(data$chr, data$pos, data$p.value, sig.level = 0.05/3403, annotate=list(ann,"tmp1"=list(col=cols[6],pch=16,cex=1.5,label=F),"tmp2"=list(col=cols[5],pch=17,cex=1,label=F)),xlab=list(label="Waist circumference",cex=2.5),ylab=list(label=expression(-log[10](italic(p))),cex=2.5))

infile <- paste0("GWAS_imputed.chr",CHR,".X49.0.0.QUALgt0.3.txt")
data <- read.table(infile, header=T, as.is=T)
data <- subset(data, AF >= maf & AF <= 1 - maf)
data <- subset(data, MISSRATE < missrate)
data <- subset(data, QUAL > qual) # imputed
data <- subset(data, !is.na(PVAL))
data <- data[,c("CHR", "cM", "PVAL")]
colnames(data) <- c("chr", "pos", "p.value")
data$source <- 1

infile <- "3cm/1/seed12345_X49.0.0.txt"
data2 <- read.table(infile, header=T)
data2 <- subset(data2, chr == CHR)
data2$pos <- (data2$start+data2$end)/2
data2 <- data2[, c("chr", "pos", "p.value")]
data2$source <- 2

infile <- paste0("3cm/1/conditional_imputed_GDS_pheno49_maf1e-04_missrate0.05_qual0.3_flank3cM_alpha",alpha,".out")
data3 <- read.table(infile, header=T, as.is=T, sep="\t")
data3 <- subset(data3, chr == CHR)
data3$pos <- (data3$start+data3$end)/2
data3 <- data3[, c("chr", "pos", "p.value.conditional")]
colnames(data3)[3] <- "p.value"
if(nrow(data3)>0) data3$source <- 3

data <- rbind(data, data2, data3)
data$p.value[data$p.value < pval.cutoff1] <- pval.cutoff1
data <- data[order(data$chr, data$pos),]
ann <- factor(data$source, levels = 1:3, labels = c("", "tmp1", "tmp2"))
p2 <- manhattan.plot(data$chr, data$pos, data$p.value, sig.level = 0.05/3403, annotate=list(ann,"tmp1"=list(col=cols[6],pch=16,cex=1.5,label=F),"tmp2"=list(col=cols[5],pch=17,cex=1,label=F)),xlab=list(label="Hip circumference",cex=2.5),ylab=list(label=expression(-log[10](italic(p))),cex=2.5))

infile <- paste0("GWAS_imputed.chr",CHR,".X21001.0.0.QUALgt0.3.txt")
data <- read.table(infile, header=T, as.is=T)
data <- subset(data, AF >= maf & AF <= 1 - maf)
data <- subset(data, MISSRATE < missrate)
data <- subset(data, QUAL > qual) # imputed
data <- subset(data, !is.na(PVAL))
data <- data[,c("CHR", "cM", "PVAL")]
colnames(data) <- c("chr", "pos", "p.value")
data$source <- 1

infile <- "3cm/1/seed12345_X21001.0.0.txt"
data2 <- read.table(infile, header=T)
data2 <- subset(data2, chr == CHR)
data2$pos <- (data2$start+data2$end)/2
data2 <- data2[, c("chr", "pos", "p.value")]
data2$source <- 2

infile <- paste0("3cm/1/conditional_imputed_GDS_pheno21001_maf1e-04_missrate0.05_qual0.3_flank3cM_alpha",alpha,".out")
data3 <- read.table(infile, header=T, as.is=T, sep="\t")
data3 <- subset(data3, chr == CHR)
data3$pos <- (data3$start+data3$end)/2
data3 <- data3[, c("chr", "pos", "p.value.conditional")]
colnames(data3)[3] <- "p.value"
if(nrow(data3)>0) data3$source <- 3

data <- rbind(data, data2, data3)
data$p.value[data$p.value < pval.cutoff1] <- pval.cutoff1
data <- data[order(data$chr, data$pos),]
ann <- factor(data$source, levels = 1:3, labels = c("", "tmp1", "tmp2"))
p5 <- manhattan.plot(data$chr, data$pos, data$p.value, sig.level = 0.05/3403, annotate=list(ann,"tmp1"=list(col=cols[6],pch=16,cex=1.5,label=F),"tmp2"=list(col=cols[5],pch=17,cex=1,label=F)),xlab=list(label="Body mass index",cex=2.5),ylab=list(label=expression(-log[10](italic(p))),cex=2.5))

infile <- paste0("GWAS_imputed.chr",CHR,".X21002.0.0.QUALgt0.3.txt")
data <- read.table(infile, header=T, as.is=T)
data <- subset(data, AF >= maf & AF <= 1 - maf)
data <- subset(data, MISSRATE < missrate)
data <- subset(data, QUAL > qual) # imputed
data <- subset(data, !is.na(PVAL))
data <- data[,c("CHR", "cM", "PVAL")]
colnames(data) <- c("chr", "pos", "p.value")
data$source <- 1

infile <- "3cm/1/seed12345_X21002.0.0.txt"
data2 <- read.table(infile, header=T)
data2 <- subset(data2, chr == CHR)
data2$pos <- (data2$start+data2$end)/2
data2 <- data2[, c("chr", "pos", "p.value")]
data2$source <- 2

infile <- paste0("3cm/1/conditional_imputed_GDS_pheno21002_maf1e-04_missrate0.05_qual0.3_flank3cM_alpha",alpha,".out")
data3 <- read.table(infile, header=T, as.is=T, sep="\t")
data3 <- subset(data3, chr == CHR)
data3$pos <- (data3$start+data3$end)/2
data3 <- data3[, c("chr", "pos", "p.value.conditional")]
colnames(data3)[3] <- "p.value"
if(nrow(data3)>0) data3$source <- 3

data <- rbind(data, data2, data3)
data$p.value[data$p.value < pval.cutoff1] <- pval.cutoff1
data <- data[order(data$chr, data$pos),]
ann <- factor(data$source, levels = 1:3, labels = c("", "tmp1", "tmp2"))
p6 <- manhattan.plot(data$chr, data$pos, data$p.value, sig.level = 0.05/3403, annotate=list(ann,"tmp1"=list(col=cols[6],pch=16,cex=1.5,label=F),"tmp2"=list(col=cols[5],pch=17,cex=1,label=F)),xlab=list(label="Body weight",cex=2.5),ylab=list(label=expression(-log[10](italic(p))),cex=2.5))

plot_grid(p1, p2, p5, p6, ncol = 1, align = "hv", labels = "AUTO", label_size=30)
ggsave("S3_Fig.tiff", width=12, height=16, bg="white")

