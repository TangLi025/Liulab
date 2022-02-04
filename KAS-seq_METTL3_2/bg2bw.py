GROUP=["KAS-seq_ALKBH3","KAS-seq_ALKBH5"]
SAMPLE=["Ctrl","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"

REFSEQ="/disk1/home/user_09/reference/annotation/mm19/mm19_Refseq.bed"

BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

rule all:
  input:
    expand("{group}/05_bedtools/bigWig/{group}_{sample}_{treatment}_{rep}_ext.bw",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("{group}/07_deeptools/computeMatrix/{group}.mat.gz",group=GROUP),
    expand("{group}/07_deeptools/plotProfile/{group}.png",group=GROUP),
    expand("{group}/07_deeptools/plotHeatmap/{group}.png",group=GROUP),
    expand("{group}/07_deeptools/plotPCA/{group}.png",group=GROUP),
    expand("{group}/07_deeptools/plotProfile/{group}_{type}.png",group=GROUP,type=["TSS","TES"]),
    
rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule bedGraphToBigWig:
  input:
    "{group}/05_bedtools/bedGraph/{group}_{sample}_{treatment}_{rep}_ext.nor.bg",
    GENOME_INDEX
  output:
    "{group}/05_bedtools/bigWig/{group}_{sample}_{treatment}_{rep}_ext.bw"
  log:
    "{group}/logs/bedGraphToBigWig/{group}_{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedGraphToBigWig {input[0]}\
      {input[1]} {output} > {log} 2>&1"

rule computeMatrix_distribution:
  input:
    bw=expand("{group}/05_bedtools/bigWig/{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix/{group}.mat.gz",
    tab="{group}/07_deeptools/computeMatrix/{group}.tab",
    bed="{group}/07_deeptools/computeMatrix/{group}.bed"
  log:
    "{group}/logs/computeMatrix_distribution/{group}.log"
  params:
    region=REFSEQ,
    blacklist=BLACKLIST,
  threads: 25
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
    mat="{group}/07_deeptools/computeMatrix/{group}.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile/{group}.png"
  log:
    "{group}/logs/plotProfile/{group}.log"
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
    mat="{group}/07_deeptools/computeMatrix/{group}.mat.gz"
  output:
    png="{group}/07_deeptools/plotHeatmap/{group}.png"
  log:
    "{group}/logs/plotHeatmap/{group}.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotHeatmap \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
      --colorMap 'Blues' \
       > {log} 2>&1"

rule computeMatrix_distribution_TSS_TES:
  input:
    bw=expand("{group}/05_bedtools/bigWig/{group}_{sample}_IP_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix/{group}_{type}.mat.gz",
    tab="{group}/07_deeptools/computeMatrix/{group}_{type}.tab",
    bed="{group}/07_deeptools/computeMatrix/{group}_{type}.bed"
  log:
    "{group}/logs/computeMatrix_distribution/{group}_{type}.log"
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
    mat="{group}/07_deeptools/computeMatrix/{group}_{type}.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile/{group}_{type}.png"
  log:
    "{group}/logs/plotProfile/{group}_{type}.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"   

rule multiBigwigSummary:
  input:
    expand("{group}/05_bedtools/bigWig/{group}_{sample}_{treatment}_{rep}_ext.bw",group=r"{group}",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    npz="{group}/07_deeptools/multiBigwigSummary/{group}.npz",
    tab="{group}/07_deeptools/multiBigwigSummary/{group}.tab"
  params:
    blacklist=BLACKLIST
  log:
    "{group}/logs/multiBigwigSummary/{group}.log"
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
    "{group}/07_deeptools/multiBigwigSummary/{group}.npz"
  output:
    "{group}/07_deeptools/plotPCA/{group}.png"
  log:
    "{group}/logs/plotPCA/{group}.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotPCA \
      -in {input} \
      --plotHeight 15 --plotWidth 10 \
      --plotFile {output} \
      > {log} 2>&1"
