SAMPLE=["DMSO","DRB","TRIP"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCh38.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh38.p13.genome.fa.fai"

rule all:
  input:
    expand("05_bedtools/bedGraph/{sample}_{treatment}_{rep}_ext.bg",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  
rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule bam2bed:
  input:
    "04_bam_unique/{sample}_{treatment}_{rep}.bam"
  output:
    "05_bedtools/bed/{sample}_{treatment}_{rep}.bed"
  log:
    "logs/bam2bed/{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools bamtobed -i {input} \
      | sort -k 1,1 1> {output} 2> {log}"
      
rule bed_extend:
  input:
    "05_bedtools/bed/{sample}_{treatment}_{rep}.bed"
  output:
    "05_bedtools/bed_extend/{sample}_{treatment}_{rep}_ext.bed"
  log:
    "logs/bed_extend/{sample}_{treatment}_{rep}.log"
  shell:
    r"""
    awk '$3-150>0 {{if ($6~"+") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$2+150,$4,$5,$6); else if ($6~"-") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$3-150,$3,$4,$5,$6)}}' \
      {input} 1> {output} 2> {log}
    """

rule bed2bedGraph:
  input:
    "05_bedtools/bed_extend/{sample}_{treatment}_{rep}_ext.bed",
    GENOME_INDEX
  output:
    "05_bedtools/bedGraph/{sample}_{treatment}_{rep}_ext.bg"
  log:
    "logs/bed2bedGraph/{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools genomecov -i {input[0]} \
      -bg -g {input[1]} \
      | /disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools sort -i 1> {output} 2> {log}"
