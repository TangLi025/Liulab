GROUP=["METTL3_3","METTL3_2"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"

REPEATS="/disk1/home/user_09/reference/annotation/mm19/mm19_repeats_family.bed"

BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

rule all:
  input:
    #expand("{group}/07_deeptools/computeMatrix_repeats/{group}.mat.gz",group=GROUP),
    #expand("{group}/07_deeptools/plotProfile_repeats/{group}.png",group=GROUP),
    #expand("{group}/07_deeptools/plotHeatmap_repeats/{group}.png",group=GROUP),
    #expand("{group}/07_deeptools/plotPCA_repeats/{group}.png",group=GROUP),
    #expand("{group}/07_deeptools/plotProfile_repeats/{group}_{type}.png",group=GROUP,type=["TSS","TES"]),
    expand("{group}/07_deeptools/plotProfile_repeats/{group}.png",group=GROUP),
    
rule computeMatrix_distribution:
  input:
    bw=expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix_repeats/{group}.mat.gz",
    tab="{group}/07_deeptools/computeMatrix_repeats/{group}.tab",
    bed="{group}/07_deeptools/computeMatrix_repeats/{group}.bed"
  log:
    "{group}/logs/computeMatrix_distribution_repeats/{group}.log"
  params:
    region=REPEATS,
    blacklist=BLACKLIST,
  threads: 100
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
    mat="{group}/07_deeptools/computeMatrix_repeats/{group}.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile_repeats/{group}.png"
  log:
    "{group}/logs/plotProfile_repeats/{group}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"
