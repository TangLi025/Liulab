TREATMENT=["input","ip"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]

SAMPLE=["p0","p10"]

GENOME_2bit="/disk/user_09/reference/genome/mm/GRCm39.genome.2bit"
GENOME="/disk/user_09/reference/genome/mm/GRCm39.genome.fa"
GENOME_INDEX="/disk/user_09/reference/genome/mm/GRCm39.genome.fa.fai"
GTF="/disk/user_09/reference/annotation/mm19/gencode.vM28.annotation.gtf"
GTF_protein_coding="/disk/user_09/reference/annotation/mm19/gencode.vM28.annotation.protein_coding.chr.gtf"
#BLACKLIST="/disk/user_09/reference/annotation/hg19/hg19_blacklist.bed"

rule all:
  input:
    "07_deeptools/plotProfile/coverage_p0p10.pdf",
    "07_deeptools/plotProfile/compare_p0p10.pdf"
      
rule faToTwoBit:
  input:
    GENOME
  output:
    GENOME_2bit
  shell:
    "/disk/user_09/anaconda3/bin/faToTwoBit {input} {output}"

rule computeMatrix_coverage:
  input:
    bw=expand("07_deeptools/bamCoverage/{sample}_{treatment}_{rep}.bw",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="07_deeptools/computeMatrix/coverage_p0p10.mat.gz",
    tab="07_deeptools/computeMatrix/coverage_p0p10.tab",
    bed="07_deeptools/computeMatrix/coverage_p0p10.bed"
  log:
    "logs/computeMatrix_distribution/coverage_p0p10.log"
  params:
    short="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p10_short.bed",
    long="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p10_long.bed",
    nc="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p10_nc.bed",
  threads: 30
  shell:
    "/disk/user_09/anaconda3/envs/deeptools/bin/computeMatrix reference-point \
      --regionsFileName {params.short} {params.long} {params.nc} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --smartLabels \
      --referencePoint center \
      --beforeRegionStartLength 3000  \
      --afterRegionStartLength 3000 \
      --binSize 10 \
      --skipZeros \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_coverage:
  input:
    mat="07_deeptools/computeMatrix/coverage_p0p10.mat.gz"
  output:
    pdf="07_deeptools/plotProfile/coverage_p0p10.pdf"
  log:
    "logs/plotProfile/coverage_p0p10.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.pdf} \
      --dpi 300 \
       > {log} 2>&1"

rule computeMatrix_compare:
  input:
    bw=expand("07_deeptools/bigwigCompare/{sample}_{rep}.bw",sample=SAMPLE,rep=REP)
  output:
    mat="07_deeptools/computeMatrix/compare_p0p10.mat.gz",
    tab="07_deeptools/computeMatrix/compare_p0p10.tab",
    bed="07_deeptools/computeMatrix/compare_p0p10.bed"
  log:
    "logs/computeMatrix_distribution/compare_p0p10.log"
  params:
    short="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p10_short.bed",
    long="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p10_long.bed",
    nc="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p10_nc.bed",
  threads: 30
  shell:
    "/disk/user_09/anaconda3/envs/deeptools/bin/computeMatrix reference-point \
      --regionsFileName {params.short} {params.long} {params.nc} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --smartLabels \
      --referencePoint center \
      --beforeRegionStartLength 3000  \
      --afterRegionStartLength 3000 \
      --binSize 10 \
      --skipZeros \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_compare:
  input:
    mat="07_deeptools/computeMatrix/compare_p0p10.mat.gz"
  output:
    pdf="07_deeptools/plotProfile/compare_p0p10.pdf"
  log:
    "logs/plotProfile/compare_p0p10.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.pdf} \
      --dpi 300 \
       > {log} 2>&1"