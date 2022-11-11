library(data.table)
p0_pAPA_closest3 <- fread("/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/02_closest/p0_pAPA_closest3.bed")
colnames(p0_pAPA_closest3) <- c("chr1","start1","end1","name1","score1","strand1",
                                "chr2","start2","end2","name2","score2","strand2",
                                "distance")

p0_pAPA_200 <- p0_pAPA_closest3[abs(p0_pAPA_closest3$distance)<=200,]
p0_pAPA_500 <- p0_pAPA_closest3[p0_pAPA_closest3$distance>-500 & p0_pAPA_closest3$distance<=0,]

p0_stream <- table(ifelse(p0_pAPA_200$distance>0,"downstream","upstream"))

p5_pAPA_closest3 <- fread("/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/02_closest/p5_pAPA_closest3.bed")
colnames(p5_pAPA_closest3) <- c("chr1","start1","end1","name1","score1","strand1",
                                 "chr2","start2","end2","name2","score2","strand2",
                                 "distance")

p5_pAPA_200 <- p5_pAPA_closest3[abs(p5_pAPA_closest3$distance)<=200,]
p5_pAPA_500 <- p5_pAPA_closest3[p5_pAPA_closest3$distance>-500 & p5_pAPA_closest3$distance<=0,]

p5_stream <- table(ifelse(p5_pAPA_200$distance>0,"downstream","upstream"))

p10_pAPA_closest3 <- fread("/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/02_closest/p10_pAPA_closest3.bed")
colnames(p10_pAPA_closest3) <- c("chr1","start1","end1","name1","score1","strand1",
                                "chr2","start2","end2","name2","score2","strand2",
                                "distance")

p10_pAPA_200 <- p10_pAPA_closest3[abs(p10_pAPA_closest3$distance)<=200,]
p10_pAPA_500 <- p10_pAPA_closest3[p10_pAPA_closest3$distance>-500 & p10_pAPA_closest3$distance<=0,]

p10_stream <- table(ifelse(p10_pAPA_200$distance>0,"downstream","upstream"))

rp2_pAPA_closest3 <- fread("/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/02_closest/rp2_pAPA_closest3.bed")
colnames(rp2_pAPA_closest3) <- c("chr1","start1","end1","name1","score1","strand1",
                                 "chr2","start2","end2","name2","score2","strand2",
                                 "distance")

rp2_pAPA_200 <- rp2_pAPA_closest3[abs(rp2_pAPA_closest3$distance)<=200,]
rp2_pAPA_500 <- rp2_pAPA_closest3[rp2_pAPA_closest3$distance>-500 & rp2_pAPA_closest3$distance<=0,]

rp2_stream <- table(ifelse(rp2_pAPA_200$distance>0,"downstream","upstream"))

fisher.test(data.frame(p0_stream,p10_stream)[,c(2,4)])
fisher.test(data.frame(p0_stream,p5_stream)[,c(2,4)])
fisher.test(data.frame(p0_stream,rp2_stream)[,c(2,4)])
fisher.test(data.frame(p5_stream,p10_stream)[,c(2,4)])
fisher.test(data.frame(p5_stream,rp2_stream)[,c(2,4)])
fisher.test(data.frame(p10_stream,rp2_stream)[,c(2,4)])

p0_pAPA_mean_distance <- c(mean(p0_pAPA_200[p0_pAPA_200$distance<0, "distance"]),
  mean(p0_pAPA_200[p0_pAPA_200$distance>0, "distance"]))
p0_pAPA_mean_distance

p5_pAPA_mean_distance <- c(mean(p5_pAPA_200[p5_pAPA_200$distance<0, "distance"]),
  mean(p5_pAPA_200[p5_pAPA_200$distance>0, "distance"]))
p5_pAPA_mean_distance

p10_pAPA_mean_distance <- c(mean(p10_pAPA_200[p10_pAPA_200$distance<0, "distance"]),
  mean(p10_pAPA_200[p10_pAPA_200$distance>0, "distance"]))
p10_pAPA_mean_distance

rp2_pAPA_mean_distance <- c(mean(rp2_pAPA_200[rp2_pAPA_200$distance<0, "distance"]),
  mean(rp2_pAPA_200[rp2_pAPA_200$distance>0, "distance"]))
rp2_pAPA_mean_distance

p0_pAPA_mean_distance <- c(mean(p0_pAPA_closest3[p0_pAPA_closest3$distance<0, "distance"]),
                           mean(p0_pAPA_closest3[p0_pAPA_closest3$distance>0, "distance"]))
p0_pAPA_mean_distance

p5_pAPA_mean_distance <- c(mean(p5_pAPA_closest3[p5_pAPA_closest3$distance<0, "distance"]),
                           mean(p5_pAPA_closest3[p5_pAPA_closest3$distance>0, "distance"]))
p5_pAPA_mean_distance

p10_pAPA_mean_distance <- c(mean(p10_pAPA_closest3[p10_pAPA_closest3$distance<0, "distance"]),
                            mean(p10_pAPA_closest3[p10_pAPA_closest3$distance>0, "distance"]))
p10_pAPA_mean_distance

rp2_pAPA_mean_distance <- c(mean(rp2_pAPA_closest3[rp2_pAPA_closest3$distance<0, "distance"]),
                            mean(rp2_pAPA_closest3[rp2_pAPA_closest3$distance>0, "distance"]))
rp2_pAPA_mean_distance

