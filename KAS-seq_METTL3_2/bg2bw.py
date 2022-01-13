GROUP=["METTL3_2"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"

REFSEQ="/disk1/home/user_09/reference/annotation/mm19/mm19_Refseq.bed"

BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

rule all:
  input:
    expand("05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
    "07_deeptools/computeMatrix/METTL3_2.mat.gz"
    
rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule bedGraphToBigWig:
  input:
    "05_bedtools/bedGraph/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.nor.bg",
    GENOME_INDEX
  output:
    "05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw"
  log:
    "logs/bedGraphToBigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedGraphToBigWig {input[0]}\
      {input[1]} {output} > {log} 2>&1"

rule computeMatrix_distribution:
  input:
    bw=expand("05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="07_deeptools/computeMatrix/{group}.mat.gz",
    tab="07_deeptools/computeMatrix/{group}.tab",
    bed="07_deeptools/computeMatrix/{group}.bed"
  log:
    "logs/computeMatrix_distribution/{group}.log"
  params:
    region=REFSEQ,
    blacklist=BLACKLIST,
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.region} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --regionBodyLength 6000 \
      --beforeRegionStartLength 3000  \
      --afterRegionStartLength 3000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --smartLabels \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"
