GROUP=["METTL3_2"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

rule all:
  input:
    expand("06_macs2/regular/KAS-seq_{group}_{sample}_{rep}_peaks.narrowPeak",group=GROUP,sample=SAMPLE,rep=REP),
    expand("06_macs2/broad/KAS-seq_{group}_{sample}_{rep}_peaks.broadPeak",group=GROUP,sample=SAMPLE,rep=REP),
    expand("06_macs2/broad_merge/KAS-seq_{group}_{sample}_peaks.broadPeak",group=GROUP,sample=SAMPLE,rep=REP)
   

rule macs2_callpeak_regular:
  input:
    "05_bedtools/bed_extend/KAS-seq_{group}_{sample}_IP_{rep}_ext.bed",
    "05_bedtools/bed_extend/KAS-seq_{group}_{sample}_input_{rep}_ext.bed"
  output:
    "06_macs2/regular/KAS-seq_{group}_{sample}_{rep}_peaks.narrowPeak",
    "06_macs2/regular/KAS-seq_{group}_{sample}_{rep}_peaks.xls"
  log:
    "logs/macs2_callpeak/regular/KAS-seq_{group}_{sample}_{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n KAS-seq_{wildcards.group}_{wildcards.sample}_{wildcards.rep} \
      -g mm -q 0.01 \
      --outdir 06_macs2/regular > {log} 2>&1"
      
rule macs2_callpeak_broad:
  input:
    "05_bedtools/bed_extend/KAS-seq_{group}_{sample}_IP_{rep}_ext.bed",
    "05_bedtools/bed_extend/KAS-seq_{group}_{sample}_input_{rep}_ext.bed"
  output:
    "06_macs2/broad/KAS-seq_{group}_{sample}_{rep}_peaks.broadPeak",
    "06_macs2/broad/KAS-seq_{group}_{sample}_{rep}_peaks.xls"
  log:
    "logs/macs2_callpeak/broad/KAS-seq_{group}_{sample}_{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n KAS-seq_{wildcards.group}_{wildcards.sample}_{wildcards.rep} \
      --broad --broad-cutoff 0.01 \
      -g mm -q 0.01 \
      --outdir 06_macs2/broad > {log} 2>&1"

rule macs2_callpeak_broad_merge:
  input:
    IP=expand("05_bedtools/bed_extend/KAS-seq_{group}_{sample}_IP_{rep}_ext.bed",group="METTL3_2",sample='{sample}',rep=REP),
    input=expand("05_bedtools/bed_extend/KAS-seq_{group}_{sample}_input_{rep}_ext.bed",group="METTL3_2",sample='{sample}',rep=REP),
  output:
    "06_macs2/broad_merge/KAS-seq_{group}_{sample}_peaks.broadPeak",
    "06_macs2/broad_merge/KAS-seq_{group}_{sample}_peaks.xls"
  log:
    "logs/macs2_callpeak/broad_merge/KAS-seq_{group}_{sample}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input.IP} -c {input.input} \
      -n KAS-seq_{wildcards.group}_{wildcards.sample} \
      --broad --broad-cutoff 0.01 \
      -g mm -q 0.01 \
      --outdir 06_macs2/broad_merge > {log} 2>&1"
