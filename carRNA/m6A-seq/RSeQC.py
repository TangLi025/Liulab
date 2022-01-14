TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]

INDEX="/disk1/home/user_09/reference/index/hisat2/hg19/hg19"
REFGENE_BED="/disk1/home/user_09/reference/annotation/hg19_GencodeCompV19.bed"


rule all:
  input:
    expand("05_FPKM_counts/HEK_{rep}.FPKM.xls",rep=REP)
    
rule FPKM_counts:
  input:
    "04_bam_raw/HEK_{rep}.bam"
  output:
    "05_FPKM_counts/HEK_{rep}.FPKM.xls"
  params:
    prefix="05_FPKM_counts/HEK_{rep}",
    refgene_bed=REFGENE_BED
  log:
    "logs/FPKM_counts/HEK_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/RSeQC/bin/FPKM_count.py \
      --input-file={input} --out-prefix={params.prefix} \
      --refgene={params.refgene_bed} 1> {log} 2>&1"
    
