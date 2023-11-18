library(SeqArray)
i <- commandArgs(T)[1]
cat("Processing chr",i,"...\n")

outdir <- "data/UKB/imputed/cm/chunk_1cM/"
setwd(outdir)

gdsfile <- paste0("data/UKB/imputed/chr",i,"_v3.impute2.gds")
mapfile <- paste0("data/UKB/imputed/cm/all.bim.gen_",i) # V1: chr; V7: cM

map <- read.table(mapfile, as.is=T)
f <- seqOpen(gdsfile)
variant.idx <- seqGetData(f, "variant.id")

for(j in seq_len(ceiling(max(map$V7)))) {
    print(j)
    idx <- which(map$V7 >= j-1 & map$V7 < j)
    if(length(idx) == 0) next
    seqSetFilter(f, variant.id = variant.idx[idx])
    seqExport(f, paste0("chr",i,".start",j-1,"cM.gds"))
}
seqClose(f)
