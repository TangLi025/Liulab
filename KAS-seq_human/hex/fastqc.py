TREATMENT=["input","IP"]
REP=["rep1","rep2"]
SAMPLE=["T00","T05","T15","T30","T60"]

rule all:
  input:
    "01_fastqc1/multiqc1/multiqc_report.html",
    "03_fastqc2/multiqc2/multiqc_report.html"
    
rule fasterq_dump:
  input:
    "00_raw_sra/{sample}_{treatment}_{rep}.sra"
  output:
    temp("00_raw_fastq/{sample}_{treatment}_{rep}.fastq")
  log:
    "logs/fasterq_dump/{sample}_{treatment}_{rep}.log"
  threads:10
  shell:
     "/disk1/home/user_09/anaconda3/envs/m6A/bin/fasterq-dump -e {threads} --split-3 -O 00_raw_fastq {input} > {log} 2>&1"

rule fastqc1:
  input:
    "00_raw_fastq/{sample}_{treatment}_{rep}.fastq"
  output:
    "01_fastqc1/{sample}_{treatment}_{rep}_fastqc.html",
    "01_fastqc1/{sample}_{treatment}_{rep}_fastqc.zip"
  log:
    "logs/fastqc1/{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o 01_fastqc1 {input} > {log} 2>&1"
    
rule multiqc1:
  input:
    expand("01_fastqc1/{sample}_{treatment}_{rep}_fastqc.html",sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("01_fastqc1/{sample}_{treatment}_{rep}_fastqc.zip",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    "01_fastqc1/multiqc1/multiqc_report.html"
  log:
    "logs/fastqc1/multiqc1.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o 01_fastqc1/multiqc1/ > {log} 2>&1"

rule trim_galore:
  input:
    "00_raw_fastq/{sample}_{treatment}_{rep}.fastq"
  output:
    "02_trim_galore/{sample}_{treatment}_{rep}_trimmed.fq"
  params:
    output_dir="02_trim_galore"
  log:
    "logs/trim_galore/{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/trim-galore/bin/trim_galore --dont_gzip \
    --length 50 -j {threads} \
    --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt \
    --output_dir {params.output_dir} {input} > {log} 2>&1"

rule fastqc2:
  input:
    "02_trim_galore/{sample}_{treatment}_{rep}_trimmed.fq"
  output:
    "03_fastqc2/{sample}_{treatment}_{rep}_trimmed_fastqc.html",
    "03_fastqc2/{sample}_{treatment}_{rep}_trimmed_fastqc.zip"
  log:
    "logs/fastqc2/{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q \
      -o 03_fastqc2 {input} > {log} 2>&1"
    
rule multiqc2:
  input:
    expand("03_fastqc2/{sample}_{treatment}_{rep}_trimmed_fastqc.html",sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("03_fastqc2/{sample}_{treatment}_{rep}_trimmed_fastqc.zip",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    "03_fastqc2/multiqc2/multiqc_report.html"
  log:
    "logs/fastqc2/multiqc2.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o 03_fastqc2/multiqc2 > {log} 2>&1"
    



