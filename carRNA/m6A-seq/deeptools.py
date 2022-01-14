TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
SAMPLE=["Lysate","Result"]
SHUFFLE=["","_shuffled"]


GENOME_2bit="/disk1/home/user_09/reference/genome/NCBIM37.genome.2bit"
GENOME="/disk1/home/user_09/reference/genome/NCBIM37.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/NCBIM37.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/mm10/gencode.vM1.annotation.gtf"

rule all:
  input:
    "07_deeptools/plotFingerprint/{sample}_m6A-seq_plotFingerprint_raw.png",
    #expand("07_deeptools/computeGCBias/{sample}_{treatment}_{rep}_raw.txt",treatment=TREATMENT,rep=REP)

rule faToTwoBit:
  input:
    GENOME
  output:
    GENOME_2bit
  shell:
    "/disk1/home/user_09/anaconda3/bin/faToTwoBit {input} {output}"

rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule plotFingerprint:
  input:
    expand("04_bam_rmdup/{sample}_{treatment}_{rep}.bam",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    "07_deeptools/plotFingerprint/{sample}_m6A-seq_plotFingerprint.png"
  log:
    "logs/plotFingerprint/{sample}_m6A-seq_plotFingerprint.log"
  params:
    labels=expand("{sample}_{treatment}_{rep}",sample=SAMPLE,treatment=["input","m6A-seq"],rep=REP),
  threads: 20
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotFingerprint \
      -b {input} -o {output} \
      --extendReads --labels {params.labels} \
      --binSize 500 \
      --numberOfProcessors {threads} > {log} 2>&1"
 
rule bamCoverage:
  input:
    bam="04_bam_rmdup/{sample}_{treatment}_{rep}.bam"
  output:
    bw="07_deeptools/bamCoverage/{sample}_{treatment}_{rep}.bw"
  log:
    "logs/bamCoverage/{sample}_{treatment}_{rep}.log"
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
      --effectiveGenomeSize 2489384235 \
      --normalizeUsing None \
      --extendReads \
      > {log} 2>&1"
      
rule computeGCBias:
  input:
    "04_bam_raw/{sample}_{treatment}_{rep}.bam"
  output:
    freq="07_deeptools/computeGCBias/{sample}_{treatment}_{rep}_raw.txt",
    png="07_deeptools/computeGCBias/{sample}_{treatment}_{rep}_raw.png"
  log:
    "logs/computeGCBias/{sample}_{treatment}_{rep}.log"
  params:
    genome=GENOME_2bit
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeGCBias \
      -b {input} \
      --effectiveGenomeSize 2489384235 \
      --genome {params.genome} \
      --GCbiasFrequenciesFile {output.freq} \
      --fragmentLength 150 \
      --biasPlot {output.png} --plotFileFormat png \
      --regionSize 300 \
      --numberOfProcessors {threads} > {log} 2>&1"
