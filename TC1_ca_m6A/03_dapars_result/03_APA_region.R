library(data.table)
PDUI_dAPA_peak_p0 <- fread("/disk/user_09/Data/03_TC1_caRNA/12_dapars/01_PDUI_dAPA_m6A_peak/PDUI_dAPA_peak_p0.bed")

sum(PDUI_dAPA_peak_p0$V7)

PDUI_dAPA_peak_p5 <- fread("/disk/user_09/Data/03_TC1_caRNA/12_dapars/01_PDUI_dAPA_m6A_peak/PDUI_dAPA_peak_p5.bed")

sum(PDUI_dAPA_peak_p5$V7)
PDUI_dAPA_peak_p10 <- fread("/disk/user_09/Data/03_TC1_caRNA/12_dapars/01_PDUI_dAPA_m6A_peak/PDUI_dAPA_peak_p10.bed")

sum(PDUI_dAPA_peak_p10$V7)
PDUI_dAPA_peak_rp2 <- fread("/disk/user_09/Data/03_TC1_caRNA/12_dapars/01_PDUI_dAPA_m6A_peak/PDUI_dAPA_peak_rp2.bed")

sum(PDUI_dAPA_peak_rp2$V7)

PDUI_dAPA_peak <- data.frame(row.names = PDUI_dAPA_peak_p0$V4,
                             p0=PDUI_dAPA_peak_p0$V7,
                             p5=PDUI_dAPA_peak_p5$V7,
                             p10=PDUI_dAPA_peak_p10$V7,
                             rp2=PDUI_dAPA_peak_rp2$V7)

################## p0p10 ######
PDUI_p0p10 <- as.data.frame(PDUI_p0p10)
rownames(PDUI_p0p10) <- substr(PDUI_p0p10$Gene,1,18)

PDUI_p0p10$dAPA_m6A_peak_gain <- PDUI_dAPA_peak[rownames(PDUI_p0p10),"p10"]-PDUI_dAPA_peak[rownames(PDUI_p0p10),"p0"]

sum(PDUI_p0p10[PDUI_p0p10$filter=="NC","dAPA_m6A_peak_gain"])
sum(PDUI_p0p10[PDUI_p0p10$filter=="UP","dAPA_m6A_peak_gain"])
sum(PDUI_p0p10[PDUI_p0p10$filter=="DOWN","dAPA_m6A_peak_gain"])
table(PDUI_p0p10[PDUI_p0p10$filter=="NC","dAPA_m6A_peak_gain"])
table(PDUI_p0p10[PDUI_p0p10$filter=="UP","dAPA_m6A_peak_gain"])
table(PDUI_p0p10[PDUI_p0p10$filter=="DOWN","dAPA_m6A_peak_gain"])

################## p0p5 ######
PDUI_p0p5 <- as.data.frame(PDUI_p0p5)
rownames(PDUI_p0p5) <- substr(PDUI_p0p5$Gene,1,18)

PDUI_p0p5$dAPA_m6A_peak_gain <- PDUI_dAPA_peak[rownames(PDUI_p0p5),"p5"]-PDUI_dAPA_peak[rownames(PDUI_p0p5),"p0"]

sum(PDUI_p0p5[PDUI_p0p5$filter=="NC","dAPA_m6A_peak_gain"])
sum(PDUI_p0p5[PDUI_p0p5$filter=="UP","dAPA_m6A_peak_gain"])
sum(PDUI_p0p5[PDUI_p0p5$filter=="DOWN","dAPA_m6A_peak_gain"])
table(PDUI_p0p5[PDUI_p0p5$filter=="NC","dAPA_m6A_peak_gain"])
table(PDUI_p0p5[PDUI_p0p5$filter=="UP","dAPA_m6A_peak_gain"])
table(PDUI_p0p5[PDUI_p0p5$filter=="DOWN","dAPA_m6A_peak_gain"])

################## p0rp2 ######
PDUI_p0rp2 <- as.data.frame(PDUI_p0rp2)
rownames(PDUI_p0rp2) <- substr(PDUI_p0rp2$Gene,1,18)

PDUI_p0rp2$dAPA_m6A_peak_gain <- PDUI_dAPA_peak[rownames(PDUI_p0rp2),"rp2"]-PDUI_dAPA_peak[rownames(PDUI_p0rp2),"p0"]

sum(PDUI_p0rp2[PDUI_p0rp2$filter=="NC","dAPA_m6A_peak_gain"])
sum(PDUI_p0rp2[PDUI_p0rp2$filter=="UP","dAPA_m6A_peak_gain"])
sum(PDUI_p0rp2[PDUI_p0rp2$filter=="DOWN","dAPA_m6A_peak_gain"])
table(PDUI_p0rp2[PDUI_p0rp2$filter=="NC","dAPA_m6A_peak_gain"])
table(PDUI_p0rp2[PDUI_p0rp2$filter=="UP","dAPA_m6A_peak_gain"])
table(PDUI_p0rp2[PDUI_p0rp2$filter=="DOWN","dAPA_m6A_peak_gain"])

################## p5p10 ######
PDUI_p5p10 <- as.data.frame(PDUI_p5p10)
rownames(PDUI_p5p10) <- substr(PDUI_p5p10$Gene,1,18)

PDUI_p5p10$dAPA_m6A_peak_gain <- PDUI_dAPA_peak[rownames(PDUI_p5p10),"p10"]-PDUI_dAPA_peak[rownames(PDUI_p5p10),"p5"]

sum(PDUI_p5p10[PDUI_p5p10$filter=="NC","dAPA_m6A_peak_gain"])
sum(PDUI_p5p10[PDUI_p5p10$filter=="UP","dAPA_m6A_peak_gain"])
sum(PDUI_p5p10[PDUI_p5p10$filter=="DOWN","dAPA_m6A_peak_gain"])
table(PDUI_p5p10[PDUI_p5p10$filter=="NC","dAPA_m6A_peak_gain"])
table(PDUI_p5p10[PDUI_p5p10$filter=="UP","dAPA_m6A_peak_gain"])
table(PDUI_p5p10[PDUI_p5p10$filter=="DOWN","dAPA_m6A_peak_gain"])


################## p5rp2 ######
PDUI_p5rp2 <- as.data.frame(PDUI_p5rp2)
rownames(PDUI_p5rp2) <- substr(PDUI_p5rp2$Gene,1,18)

PDUI_p5rp2$dAPA_m6A_peak_gain <- PDUI_dAPA_peak[rownames(PDUI_p5rp2),"rp2"]-PDUI_dAPA_peak[rownames(PDUI_p5rp2),"p5"]

sum(PDUI_p5rp2[PDUI_p5rp2$filter=="NC","dAPA_m6A_peak_gain"])
sum(PDUI_p5rp2[PDUI_p5rp2$filter=="UP","dAPA_m6A_peak_gain"])
sum(PDUI_p5rp2[PDUI_p5rp2$filter=="DOWN","dAPA_m6A_peak_gain"])
table(PDUI_p5rp2[PDUI_p5rp2$filter=="NC","dAPA_m6A_peak_gain"])
table(PDUI_p5rp2[PDUI_p5rp2$filter=="UP","dAPA_m6A_peak_gain"])
table(PDUI_p5rp2[PDUI_p5rp2$filter=="DOWN","dAPA_m6A_peak_gain"])


################## p10rp2 ######
PDUI_p10rp2 <- as.data.frame(PDUI_p10rp2)
rownames(PDUI_p10rp2) <- substr(PDUI_p10rp2$Gene,1,18)

PDUI_p10rp2$dAPA_m6A_peak_gain <- PDUI_dAPA_peak[rownames(PDUI_p10rp2),"rp2"]-PDUI_dAPA_peak[rownames(PDUI_p10rp2),"p10"]

sum(PDUI_p10rp2[PDUI_p10rp2$filter=="NC","dAPA_m6A_peak_gain"])
sum(PDUI_p10rp2[PDUI_p10rp2$filter=="UP","dAPA_m6A_peak_gain"])
sum(PDUI_p10rp2[PDUI_p10rp2$filter=="DOWN","dAPA_m6A_peak_gain"])
table(PDUI_p10rp2[PDUI_p10rp2$filter=="NC","dAPA_m6A_peak_gain"])
table(PDUI_p10rp2[PDUI_p10rp2$filter=="UP","dAPA_m6A_peak_gain"])
table(PDUI_p10rp2[PDUI_p10rp2$filter=="DOWN","dAPA_m6A_peak_gain"])
