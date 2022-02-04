TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]

GENOME_2bit="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.2bit"
GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/gencode.v19.annotation.gtf"
BLACKLIST="/disk1/home/user_09/reference/annotation/hg19_blacklist.bed"

rule all:
  input:
    expand("07_deeptools/bamCoverage/HEK_{treatment}_{rep}.bw",treatment=TREATMENT,rep=REP),
    #expand("07_deeptools/plotFingerprint/HEK_KAS-seq_plotFingerprint_{dup}.png",dup=DUP),
    #expand("07_deeptools/computeGCBias/HEK_{treatment}_{rep}_{dup}.txt",treatment=TREATMENT,rep=REP,dup=DUP),
    #expand("07_deeptools/multiBamSummary/HEK_{dup}.npz",dup=DUP),
    #expand("07_deeptools/plotCorrelation/HEK_{dup}.png",dup=DUP),
    #expand("07_deeptools/plotCorrelation/HEK_{dup}.tab",dup=DUP),
    #"07_deeptools/computeMatrix/HEK.mat.gz",
    #"07_deeptools/computeMatrix/HEK.tab",
    #"07_deeptools/computeMatrix/HEK.bed",
    #"07_deeptools/plotProfile/HEK.png"
    
    
rule faToTwoBit:
  input:
    GENOME
  output:
    GENOME_2bit
  shell:
    "/disk1/home/user_09/anaconda3/bin/faToTwoBit {input} {output}"

rule bamCoverage:
  input:
    bam="04_bam_rmdup/HEK_{treatment}_{rep}.bam"
  output:
    bw="07_deeptools/bamCoverage/HEK_{treatment}_{rep}.bw"
  log:
    "bamCoverage/HEK_{treatment}_{rep}.bam"
  params:
    blacklist=BLACKLIST
  threads:40
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/bamCoverage \
      --bam {input.bam} --outFileName {output.bw} \
      --outFileFormat bigwig \
      --binSize 50 \
      --blackListFileName {params.blacklist} \
      --numberOfProcessors {threads} \
      --effectiveGenomeSize 2862010578 \
      --normalizeUsing None \
      --extendReads 150 \
      > {log} 2>&1"

rule plotFingerprint:
  input:
    "04_bam_{dup}/HEK_input_rep1.bam",
    "04_bam_{dup}/HEK_input_rep2.bam",
    "04_bam_{dup}/HEK_IP_rep1.bam",
    "04_bam_{dup}/HEK_IP_rep2.bam"
  output:
    "07_deeptools/plotFingerprint/HEK_KAS-seq_plotFingerprint_{dup}.png"
  log:
    "logs/plotFingerprint/HEK_KAS-seq_plotFingerprint_{dup}.log"
  params:
    labels=expand("{sample}_{rep}",sample=["input","KAS-seq"],rep=REP),
    title="""Fingerprints of KAS-seq data"""
  threads: 25
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotFingerprint \
      -b {input} -o {output} \
      --extendReads 150 --labels {params.labels} \
      --binSize 500 \
      --numberOfProcessors {threads} > {log} 2>&1"
      
rule computeGCBias:
  input:
    bam="04_bam_{dup}/HEK_{treatment}_{rep}.bam",
    genome=GENOME_2bit
  output:
    freq="07_deeptools/computeGCBias/HEK_{treatment}_{rep}_{dup}.txt",
    png="07_deeptools/computeGCBias/HEK_{treatment}_{rep}_{dup}.png"
  log:
    "logs/computeGCBias/HEK_{treatment}_{rep}_{dup}.log"
  params:
    genome=GENOME_2bit
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeGCBias \
      -b {input.bam} \
      --effectiveGenomeSize 2862010578 \
      --genome {input.genome} \
      --GCbiasFrequenciesFile {output.freq} \
      --fragmentLength 150 \
      --biasPlot {output.png} --plotFileFormat png \
      --regionSize 300 \
      --numberOfProcessors {threads} > {log} 2>&1"
      
rule multiBamSummary:
  input:
    "04_bam_{dup}/HEK_IP_rep1.bam",
    "04_bam_{dup}/HEK_IP_rep2.bam"
  output:
    npz="07_deeptools/multiBamSummary/HEK_{dup}.npz"
  params:
    blacklist=BLACKLIST
  log:
    "logs/multiBamSummary/HEK_{dup}.log"
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/multiBamSummary bins \
      -b {input} \
      --outFileName {output.npz} \
      --binSize 10000 \
      --numberOfProcessors {threads} \
      --extendReads 150 \
      --blackListFileName {params.blacklist} \
       > {log} 2>&1"

rule plotCorrelation:
  input:
    npz="07_deeptools/multiBamSummary/HEK_{dup}.npz"
  output:
    png="07_deeptools/plotCorrelation/HEK_{dup}.png",
    tab="07_deeptools/plotCorrelation/HEK_{dup}.tab"
  log:
    "logs/plotCorrelation/HEK_{dup}.log"
  params:
    labels=["rep1","rep2"]
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotCorrelation \
      --corData {input.npz} \
      --corMethod pearson \
      --whatToPlot scatterplot \
      --plotFile {output.png} \
      --labels {params.labels} \
      --outFileCorMatrix {output.tab} \
      --xRange 0 100 \
      --yRange 0 100 \
      --removeOutliers \
       > {log} 2>&1"
       
rule computeMatrix_distribution:
  input:
    bw=expand("05_bedtools/HEK_{treatment}_{rep}.bw",treatment=TREATMENT,rep=REP)
  output:
    mat="07_deeptools/computeMatrix/HEK.mat.gz",
    tab="07_deeptools/computeMatrix/HEK.tab",
    bed="07_deeptools/computeMatrix/HEK.bed"
  log:
    "logs/computeMatrix_distribution/HEK.log"
  params:
    gtf=GTF,
    blacklist=BLACKLIST,
    labels=["input_rep1","input_rep2","IP_rep1","IP_rep2"]
  threads: 30
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --regionBodyLength 10000 \
      --beforeRegionStartLength 3000  \
      --afterRegionStartLength 3000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --samplesLabel {params.labels} \
      --numberOfProcessors {threads} \
       > {log} 2>&1"

rule plotProfile:
  input:
    mat="07_deeptools/computeMatrix/HEK.mat.gz"
  output:
    png="07_deeptools/plotProfile/HEK.png"
  log:
    "logs/plotProfile/HEK.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"
