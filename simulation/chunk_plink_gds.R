library(SeqArray)
gdsfile <- "data/UKB/plink/ukb_hap_merged.gds"
outdir <- "data/UKB/plink/chunk_5.1cM/"
setwd(outdir)

map <- NULL
for(i in 1:22) {
    tmp <- read.table(paste0("ukb_",i,".rMap"))
    tmp$chr <- i
    map <- rbind(map, tmp)
}

f <- seqOpen(gdsfile)
variant.idx <- seqGetData(f, "variant.id")
for(i in 1:22) {
    print(i)
    for(j in seq_len(ceiling(max(map$V2[map$chr == i])))) {
        print(j)
        idx <- which(map$V2 >= j-1 & map$V2 < j)
	if(length(idx) == 0) next
	seqSetFilter(f, variant.id = variant.idx[idx])
	seqExport(f, paste0("chr",i,".start",j-1,"cM.gds"))
    }
}
seqClose(f)
