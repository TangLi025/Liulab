SAMPLE=["Lysate","Result"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","dedup"]
STRAND=["pos","neg"]


rule all:
  input:
    expand("06_macs2_bg/{dup}/{sample}_{rep}_{strand}_peaks.xls",dup=DUP,sample=SAMPLE,rep=REP,strand=STRAND)

rule macs2_callpeak:
  input:
    "05_bam_{dup}_separated/{sample}_IP_{rep}_{strand}.bam",
    "05_bam_{dup}_separated/{sample}_input_{rep}_{strand}.bam"
  output:
    "06_macs2_bg/{dup}/{sample}_{rep}_{strand}_peaks.xls",
    "06_macs2_bg/{dup}/{sample}_{rep}_{strand}_peaks.narrowPeak"
  log:
    "logs/macs2_callpeak_bg/{dup}/{sample}_{rep}_{strand}.log"
  params:
    out_name="{sample}_{rep}_{strand}",
    out_dir="06_macs2_bg/{dup}"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n {params.out_name} \
      -f BAM --verbose 3 --nomodel --extsize 150 \
      --keep-dup 5 \
      -B --SPMR \
      -g 1.3e8 -q 0.01 \
      --outdir {params.out_dir} > {log} 2>&1"
