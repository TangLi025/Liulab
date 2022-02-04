
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

sample=c("KAS-seq","Random")

# R packages
library(tidyverse)
library(ggplot2)
library(reshape2)

# read annotation file from macs2 and homer
anno1 <- read.delim(snakemake@input[[1]],header = T)
anno2 <- read.delim(snakemake@input[[2]],header = T)

# save as list
df <- list(Lysate = anno1$Annotation,Result = anno2$Annotation)

# save separately
lapply(1:2, function(x){
  tmp = sapply(strsplit(df[[x]],split = '\\('),'[',1)
  res = table(tmp)
}) %>% do.call('rbind',.) %>% as.data.frame() -> type_comb

# add name
type_comb$name = sample

# width to length
final <- melt(type_comb)

# factor
final$name <- factor(final$name,levels = sample)

# draw pie
pdf(file=snakemake@output[[1]])
ggplot(final,aes(x = '',y = value,fill = variable)) +
  geom_col(position = position_fill()) +
  theme_void() +
  theme(legend.position = 'bottom',
        strip.text.x = element_text(size= 20)) +
  facet_wrap(~name) +
  coord_polar(theta = 'y') +
  scale_fill_brewer(palette = 'Set1',name = 'Region types')
dev.off()
# draw bar
# width to length
final <- melt(type_comb)

percent <- type_comb
rownames(percent) <- percent$name
percent <- percent[,-6]
percent[sample[1],] <- percent[sample[1],]/sum(percent[sample[1],])*100
percent[sample[2],] <- percent[sample[2],]/sum(percent[sample[2],])*100
write.table(percent,file = snakemake@output[[3]])
# draw bar
pdf(file=snakemake@output[[2]])
ggplot(final,aes(x = name,y = value,fill = variable)) +
  geom_col(position = position_fill()) +
  theme_bw(base_size = 18) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = 'right',
        strip.text.x = element_text(size= 18)) +
  scale_fill_brewer(palette = 'Set1',name = 'Region types') +
  xlab('') + ylab('Percent')
dev.off()