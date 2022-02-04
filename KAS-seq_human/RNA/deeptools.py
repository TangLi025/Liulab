TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
SHUFFLE=["","_shuffled"]


GENOME_2bit="/disk1/home/user_09/reference/genome/GRCh38.p13.genome.2bit"
GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/gencode.v19.annotation.gtf"

rule all:
  input:
    "07_deeptools/plotFingerprint/HEK_KAS-seq_plotFingerprint_raw.png",
    expand("07_deeptools/computeGCBias/HEK_{treatment}_{rep}_raw.txt",treatment=TREATMENT,rep=REP)
    GENOME_2bit

rule faToTwoBit:
  input:
    GENOME
  output:
    GENOME_2bit
  shell:
    "/disk1/home/user_09/anaconda3/bin/faToTwoBit {input} {output}"


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
