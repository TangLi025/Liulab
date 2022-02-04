TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
#SAMPLE=["DMSO"]
SAMPLE=["T00","T05","T15","T30","T60"]
MERGE=["","_merge"]

GENOME_2bit="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.2bit"
GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.gtf"
GTF_protein_coding="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.gene_coding.gtf"
BLACKLIST="/disk1/home/user_09/reference/annotation/hg19/hg19_blacklist.bed"

SHUFFLE=["","_shuffled"]

rule all:
  input:
    #"07_deeptools/plotFingerprint/HEK_KAS-seq_plotFingerprint_raw.png",
    #expand("07_deeptools/computeGCBias/HEK_{treatment}_{rep}_raw.txt",treatment=TREATMENT,rep=REP)
    #GENOME_2bit
    expand("07_deeptools/bamCoverage/{sample}_input_merge.bw",sample=SAMPLE),
    "07_deeptools/plotProfile/computeMatrix_reference_point.png"

rule faToTwoBit:
  input:
    GENOME
  output:
    GENOME_2bit
  shell:
    "/disk1/home/user_09/anaconda3/bin/faToTwoBit {input} {output}"

rule bam_merge:
  input:
    "04_bam_rmdup/{sample}_{treatment}_rep1.bam",
    "04_bam_rmdup/{sample}_{treatment}_rep2.bam"
  output:
    "04_bam_rmdup/{sample}_{treatment}_merge.bam"
  log:
    "logs/bam_merge/{sample}_{treatment}.log"
  threads:4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools merge -f -@ {threads} {output} {input[0]} {input[1]} > {log} 2>&1"

rule bam_index:
  input:
    "04_bam_rmdup/{sample}_{treatment}_merge.bam"
  output:
    "04_bam_rmdup/{sample}_{treatment}_merge.bam.bai"
  log:
    "logs/bam_index_merge/{sample}_{treatment}_merge.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"

rule bamCoverage:
  input:
    bam="04_bam_rmdup/{sample}_{treatment}_{rep}.bam",
    index="04_bam_rmdup/{sample}_{treatment}_{rep}.bam.bai"
  output:
    bw="07_deeptools/bamCoverage/{sample}_{treatment}_{rep}.bw"
  log:
    "logs/bamCoverage/{sample}_{treatment}_{rep}.log"
  params:
    blacklist=BLACKLIST
  threads:10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/bamCoverage \
      --bam {input.bam} --outFileName {output.bw} \
      --outFileFormat bigwig \
      --binSize 50 \
      --normalizeUsing None \
      --numberOfProcessors {threads} \
      --effectiveGenomeSize 2862010578 \
      --extendReads 150 \
      > {log} 2>&1"
      
rule computeMatrix_reference_point:
  input:
    bw=expand("07_deeptools/bamCoverage/{sample}_input_merge.bw",sample=SAMPLE)
  output:
    mat="07_deeptools/computeMatrix/computeMatrix_reference_point.mat.gz",
    tab="07_deeptools/computeMatrix/computeMatrix_reference_point.tab",
    bed="07_deeptools/computeMatrix/computeMatrix_reference_point.bed"
  log:
    "logs/computeMatrix_distribution/computeMatrix_reference_point_HEK.log"
  params:
    gtf=GTF_protein_coding,
    blacklist=BLACKLIST,
    labels=SAMPLE
  threads: 45
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix reference-point \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --referencePoint TSS \
      --beforeRegionStartLength 2000  \
      --afterRegionStartLength 2000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --samplesLabel {params.labels} \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_reference_point:
  input:
    mat="07_deeptools/computeMatrix/computeMatrix_reference_point.mat.gz"
  output:
    png="07_deeptools/plotProfile/computeMatrix_reference_point.png"
  log:
    "logs/plotProfile/computeMatrix_reference_point.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"


rule plotFingerprint:
  input:
    expand("04_bam_raw/HEK_{treatment}_{rep}.bam",treatment=TREATMENT,rep=REP)
  output:
    "07_deeptools/plotFingerprint/HEK_KAS-seq_plotFingerprint_raw.png"
  log:
    "logs/plotFingerprint/HEK_KAS-seq_plotFingerprint.log"
  params:
    labels=expand("{sample}_{rep}",sample=["input","KAS-seq"],rep=REP),
    title="""Fingerprints of KAS-seq data"""
  threads: 40
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotFingerprint \
      -b {input} -o {output} \
      --extendReads 150 --labels {params.labels} \
      --binSize 500 \
      --numberOfProcessors {threads} > {log} 2>&1"
      
rule computeGCBias:
  input:
    "04_bam_raw/HEK_{treatment}_{rep}.bam"
  output:
    freq="07_deeptools/computeGCBias/HEK_{treatment}_{rep}_raw.txt",
    png="07_deeptools/computeGCBias/HEK_{treatment}_{rep}_raw.png"
  log:
    "logs/computeGCBias/HEK_{treatment}_{rep}.log"
  params:
    genome=GENOME_2bit
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeGCBias \
      -b {input} \
      --effectiveGenomeSize 2862010578 \
      --genome {params.genome} \
      --GCbiasFrequenciesFile {output.freq} \
      --fragmentLength 150 \
      --biasPlot {output.png} --plotFileFormat png \
      --regionSize 300 \
      --numberOfProcessors {threads} > {log} 2>&1"
