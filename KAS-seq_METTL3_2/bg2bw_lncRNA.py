GROUP=["METTL3_1","METTL3_3","METTL14","YTHDC1"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"

REFSEQ="/disk1/home/user_09/reference/annotation/mm19/gencode.vM28.annotation.lncRNA.gtf"

BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

rule all:
  input:
    expand("{group}/07_deeptools/computeMatrix/{group}_lncRNA.mat.gz",group=GROUP),
    expand("{group}/07_deeptools/plotProfile/{group}_lncRNA.png",group=GROUP),
    expand("{group}/07_deeptools/plotHeatmap/{group}_lncRNA.png",group=GROUP),
    expand("{group}/07_deeptools/plotPCA/{group}_lncRNA.png",group=GROUP),
    expand("{group}/07_deeptools/plotProfile/{group}_{type}_lncRNA.png",group=GROUP,type=["TSS","TES"])
    
    
rule computeMatrix_distribution_TSS_TES:
  input:
    bw=expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_IP_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix/{group}_{type}_lncRNA.mat.gz",
    tab="{group}/07_deeptools/computeMatrix/{group}_{type}_lncRNA.tab",
    bed="{group}/07_deeptools/computeMatrix/{group}_{type}_lncRNA.bed"
  log:
    "{group}/logs/computeMatrix_distribution/{group}_{type}_lncRNA.log"
  params:
    type=r'{type}',
    region=REFSEQ,
    blacklist=BLACKLIST,
  threads: 20
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix reference-point \
      --regionsFileName {params.region} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --referencePoint {params.type} \
      --beforeRegionStartLength 5000  \
      --afterRegionStartLength 5000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --smartLabels \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_TSS_TES:
  input:
    mat="{group}/07_deeptools/computeMatrix/{group}_{type}_lncRNA.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile/{group}_{type}_lncRNA.png"
  log:
    "{group}/logs/plotProfile/{group}_{type}_lncRNA.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"   
    
rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule bedGraphToBigWig:
  input:
    "{group}/05_bedtools/bedGraph/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.nor.bg",
    GENOME_INDEX
  output:
    "{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw"
  log:
    "{group}/logs/bedGraphToBigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedGraphToBigWig {input[0]}\
      {input[1]} {output} > {log} 2>&1"

rule computeMatrix_distribution:
  input:
    bw=expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix/{group}_lncRNA.mat.gz",
    tab="{group}/07_deeptools/computeMatrix/{group}_lncRNA.tab",
    bed="{group}/07_deeptools/computeMatrix/{group}_lncRNA.bed"
  log:
    "{group}/logs/computeMatrix_distribution/{group}_lncRNA.log"
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

rule plotProfile:
  input:
    mat="{group}/07_deeptools/computeMatrix/{group}_lncRNA.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile/{group}_lncRNA.png"
  log:
    "{group}/logs/plotProfile/{group}_lncRNA.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"
       
rule plotHeatmap:
  input:
    mat="{group}/07_deeptools/computeMatrix/{group}_lncRNA.mat.gz"
  output:
    png="{group}/07_deeptools/plotHeatmap/{group}_lncRNA.png"
  log:
    "{group}/logs/plotHeatmap/{group}_lncRNA.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotHeatmap \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
      --colorMap 'Blues' \
       > {log} 2>&1"

rule multiBigwigSummary:
  input:
    expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=r"{group}",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    npz="{group}/07_deeptools/multiBigwigSummary/{group}_lncRNA.npz",
    tab="{group}/07_deeptools/multiBigwigSummary/{group}_lncRNA.tab"
  params:
    blacklist=BLACKLIST
  log:
    "{group}/logs/multiBigwigSummary/{group}_lncRNA.log"
  threads: 15
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/multiBigwigSummary bins \
      -b {input} \
      --outFileName {output.npz} \
      --outRawCounts {output.tab} \
      --binSize 10000 \
      --numberOfProcessors {threads} \
      --blackListFileName {params.blacklist} \
       > {log} 2>&1"

rule plotPCA:
  input:
    "{group}/07_deeptools/multiBigwigSummary/{group}_lncRNA.npz"
  output:
    "{group}/07_deeptools/plotPCA/{group}_lncRNA.png"
  log:
    "{group}/logs/plotPCA/{group}_lncRNA.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotPCA \
      -in {input} \
      --plotHeight 15 --plotWidth 10 \
      --plotFile {output} \
      > {log} 2>&1"
