GROUP=["METTL3_2"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]


rule all:
  input:
    "01_fastqc1/multiqc_report.html",
    "03_fastqc2/multiqc_report.html"
    
rule fastqc1:
  input:
    "00_raw_fastq/KAS-seq_{group}_{sample}_{treatment}_{rep}.fastq.gz"
  output:
    "01_fastqc1/KAS-seq_{group}_{sample}_{treatment}_{rep}_fastqc.html",
    "01_fastqc1/KAS-seq_{group}_{sample}_{treatment}_{rep}_fastqc.zip"
  log:
    "logs/fastqc1/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 2
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o 01_fastqc1 {input} > {log} 2>&1"
    
rule multiqc1:
  input:
    expand("01_fastqc1/KAS-seq_{group}_{sample}_{treatment}_{rep}_fastqc.html",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("01_fastqc1/KAS-seq_{group}_{sample}_{treatment}_{rep}_fastqc.zip",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    "01_fastqc1/multiqc_report.html"
  log:
    "logs/fastqc1/multiqc1.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o 01_fastqc1/ > {log} 2>&1"

rule trim_galore:
  input:
    "00_raw_fastq/KAS-seq_{group}_{sample}_{treatment}_{rep}.fastq.gz"
  output:
    "02_trim_galore/KAS-seq_{group}_{sample}_{treatment}_{rep}_trimmed.fq.gz"
  params:
    output_dir="02_trim_galore"
  log:
    "logs/trim_galore/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/trim-galore/bin/trim_galore \
    --length 50 -j {threads} \
    --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt \
    --output_dir {params.output_dir} {input} > {log} 2>&1"

rule fastqc2:
  input:
    "02_trim_galore/KAS-seq_{group}_{sample}_{treatment}_{rep}_trimmed.fq.gz"
  output:
    "03_fastqc2/KAS-seq_{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.html",
    "03_fastqc2/KAS-seq_{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.zip"
  log:
    "logs/fastqc2/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 2
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q \
      -o 03_fastqc2 {input} > {log} 2>&1"
    
rule multiqc2:
  input:
    expand("03_fastqc2/KAS-seq_{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.html",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("03_fastqc2/KAS-seq_{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.zip",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    "03_fastqc2/multiqc_report.html"
  log:
    "logs/fastqc2/multiqc2.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o 03_fastqc2/ > {log} 2>&1"
    



