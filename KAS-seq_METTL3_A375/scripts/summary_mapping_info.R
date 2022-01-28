log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")


library(tidyverse)
      
# 提取信息
lapply(snakemake@input,function(x){
  # 读取文件
  align_info <- read.delim(x,header = F,fill = T)
  # 提取唯一比对率
  unique_map <- sapply(strsplit(align_info$V1[4],split = "\\(|\\)"),'[',2)
  # 提取总比对率
  total_map <- substr(align_info$V1[15],1,6)
  # 样本名称
  name <- substr(x,1,11)
  # 合并
  res <- cbind(name,unique_map,total_map)
}) %>% Reduce('rbind',.) %>%
  as.data.frame() -> mapping_info

path <- snakemake@output
write.table(mapping_info,file = path)