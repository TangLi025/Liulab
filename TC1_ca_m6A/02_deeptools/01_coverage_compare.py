TREATMENT=["input","ip"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]

SAMPLE=["p0","p5","p10","rp2"]

GENOME_2bit="/disk/user_09/reference/genome/mm/GRCm39.genome.2bit"
GENOME="/disk/user_09/reference/genome/mm/GRCm39.genome.fa"
GENOME_INDEX="/disk/user_09/reference/genome/mm/GRCm39.genome.fa.fai"
GTF="/disk/user_09/reference/annotation/mm19/gencode.vM28.annotation.gtf"
GTF_protein_coding="/disk/user_09/reference/annotation/mm19/gencode.vM28.annotation.protein_coding.chr.gtf"
#BLACKLIST="/disk/user_09/reference/annotation/hg19/hg19_blacklist.bed"

rule all:
  input:
    expand("07_deeptools/bigwigCompare/{sample}_{rep}.bw",sample=SAMPLE,rep=REP)

    
rule faToTwoBit:
  input:
    GENOME
  output:
    GENOME_2bit
  shell:
    "/disk/user_09/anaconda3/bin/faToTwoBit {input} {output}"

rule bamCoverage:
  input:
    "/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/{sample}_{treatment}_{rep}.bam"
  output:
    bw="07_deeptools/bamCoverage/{sample}_{treatment}_{rep}.bw"
  log:
    "logs/bamCoverage/{sample}_{treatment}_{rep}.log"
  threads: 10
  shell:
    "/disk/user_09/anaconda3/envs/deeptools/bin/bamCoverage \
      -b {input} \
      --outFileName {output.bw} \
      --outFileFormat bigwig \
      --numberOfProcessors {threads} \
      --effectiveGenomeSize 2494787188 \
      --normalizeUsing RPKM \
       > {log} 2>&1"

rule bigwigCompare:
  input:
    ip="07_deeptools/bamCoverage/{sample}_ip_{rep}.bw",
    input="07_deeptools/bamCoverage/{sample}_input_{rep}.bw"
  output:
    bw="07_deeptools/bigwigCompare/{sample}_{rep}.bw"
  log:
    "logs/bigwigCompare/{sample}_{rep}.log"
  threads: 10
  shell:
    "/disk/user_09/anaconda3/envs/deeptools/bin/bigwigCompare  \
      -b1 {input.ip} \
      -b2 {input.input} \
      --skipZeroOverZero \
      --skipNAs \
      --operation log2 \
      --outFileName {output.bw} \
      --outFileFormat bigwig \
      --numberOfProcessors {threads} \
       > {log} 2>&1"
