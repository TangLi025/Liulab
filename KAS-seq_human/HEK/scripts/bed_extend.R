log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

bed <- read.table(snakemake@input[[1]])
chrom_size <- read.table(snakemake@input[[2]],row.names = 1)

for (i in c(1:dim(bed)[1])){
  if (bed[i,6]=="+"){
    if (bed[i,2]+150 <= chrom_size[bed[i,1],1]){
      bed[i,3] <- bed[i,2]+150
    }
    else{
      bed[i,3] <- chrom_size[bed[i,1],1]
    }
  }
  else {
    if (bed[i,3] >150){
      bed[i,2] <- bed[i,3] - 150
    }
    else{
      bed[i,2] <- 1
    }
  }
}

write.table(bed,snakemake@output,quote = FALSE,sep = "\t",
            row.names = FALSE,col.names = FALSE)
