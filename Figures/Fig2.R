# Fig 2. Distribution of average numbers of non-zero IBD sharing per individual in local IBD matrices.

library(ggplot2)

d3 <- read.table("3cm/1/nonzero_3cm.txt", header=T)
d3$cm <- "RaPID 3cM"
d5 <- read.table("5cm/1/nonzero_5cm.txt", header=T)
d5$cm <- "RaPID 5cM"
d10 <- read.table("10cm/1/nonzero_10cm.txt", header=T)
d10$cm <- "RaPID 10cM"
hapIBD.d3 <- read.table("hapIBD/3cm/nonzero_3cm.txt", header=T)
hapIBD.d3$cm <- "hap-IBD 3cM"
hapIBD.d5 <- read.table("hapIBD/5cm/nonzero_5cm.txt", header=T)
hapIBD.d5$cm <- "hap-IBD 5cM"
hapIBD.d10 <- read.table("hapIBD/10cm/nonzero_10cm.txt", header=T)
hapIBD.d10$cm <- "hap-IBD 10cM"
d <- rbind(d3, d5, d10, hapIBD.d3, hapIBD.d5, hapIBD.d10)
d$Length <- factor(d$cm, levels = c("RaPID 3cM", "hap-IBD 3cM", "RaPID 5cM", "hap-IBD 5cM", "RaPID 10cM", "hap-IBD 10cM"))

ggplot(d, aes(x=ratio, y=after_stat(scaled), color=Length, linetype=Length)) +
    geom_density(show.legend = FALSE) +
    stat_density(geom = "line", position = "identity") +
    scale_x_log10(breaks = c(1,2,5,10,20,50,100,200,500)) +
    labs(x = "Non-zero off-diagonal elements / N")

ggsave("Fig2.eps", width=6, height=3)
