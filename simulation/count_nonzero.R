ids <- read.table("ids_after_filtering")[,1]
ids_exclude <- read.table("withdrawal.csv")[,1]
ids <- ids[!ids %in% ids_exclude]

indir <- commandArgs(T)[1]
# 3cm/1/chunk.1cM/
# 5cm/1/chunk.1cM/
# 10cm/1/chunk.1cM/
# hapIBD/3cm/chunk.1cM/
# hapIBD/5cm/chunk.1cM/
# hapIBD/10cm/chunk.1cM/
outfile <- commandArgs(T)[2]
# 3cm/1/nonzero_3cm.txt
# 5cm/1/nonzero_5cm.txt
# 10cm/1/nonzero_10cm.txt
# hapIBD/3cm/nonzero_3cm.txt
# hapIBD/5cm/nonzero_5cm.txt
# hapIBD/10cm/nonzero_10cm.txt

out <- rbind(cbind(1,0:267),cbind(2,0:251),cbind(3,0:218),cbind(4,0:202),cbind(5,0:197),cbind(6,0:186),cbind(7,0:178),cbind(8,0:161),cbind(9,0:157),cbind(10,0:169),cbind(11,0:154),cbind(12,0:165),cbind(13,0:127),cbind(14,0:116),cbind(15,0:117),cbind(16,0:126),cbind(17,0:129),cbind(18,0:116),cbind(19,0:106),cbind(20,0:107),cbind(21,0:62),cbind(22,0:70))
out <- as.data.frame(out)
names(out) <- c("chr", "start")
out$N <- rep(NA, nrow(out))
out$non.zero <- rep(NA, nrow(out))
library(Matrix)
for(i in 1:nrow(out)) {
    ibd <- get(load(paste0(indir, "/chr",out$chr[i],".start",out$start[i],"cM.RData")))
    ibd$matrix <- ibd$matrix[rownames(ibd$matrix) %in% ids, colnames(ibd$matrix) %in% ids]
    out$N[i] <- nrow(ibd$matrix)
    out$non.zero[i] <- length(ibd$matrix@i)
}
out$ratio <- out$non.zero/out$N
write.table(out, outfile, quote=F, row.names=F,col.names=T)
