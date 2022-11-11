SAMPLE=["CTRL","SENE"]

REP=["rep1","rep2"]


rule all:
  input:
    "01_fastqc1/multiqc_report.html",
    "03_fastqc2/multiqc_report.html"
    
rule fastqc1:
  input:
    "00_raw_fastq/{sample}_{rep}.fastq"
  output:
    "01_fastqc1/{sample}_{rep}_fastqc.html",
    "01_fastqc1/{sample}_{rep}_fastqc.zip"
  log:
    "logs/fastqc1/{sample}_{rep}.log"
  params:
    out_dir="01_fastqc1"
  threads: 2
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc1:
  input:
    expand("01_fastqc1/{sample}_{rep}_fastqc.html",sample=SAMPLE,rep=REP),
    expand("01_fastqc1/{sample}_{rep}_fastqc.zip",sample=SAMPLE,rep=REP)
  output:
    "01_fastqc1/multiqc_report.html"
  log:
    "logs/fastqc1/multiqc1.log"
  threads: 1
  params:
    out_dir="01_fastqc1"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o {params.out_dir} > {log} 2>&1"

rule trim_galore:
  input:
    "00_raw_fastq/{sample}_{rep}.fastq"
  output:
    "02_trim_galore/{sample}_{rep}_trimmed.fq"
  params:
    output_dir="02_trim_galore"
  log:
    "logs/trim_galore/{sample}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/trim-galore/bin/trim_galore \
    -j {threads} \
    --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt \
    --output_dir {params.output_dir} {input} > {log} 2>&1"

rule fastqc2:
  input:
    "02_trim_galore/{sample}_{rep}_trimmed.fq"
  output:
    "03_fastqc2/{sample}_{rep}_trimmed_fastqc.html",
    "03_fastqc2/{sample}_{rep}_trimmed_fastqc.zip"
  log:
    "logs/fastqc2/{sample}_{rep}.log"
  threads: 2
  params:
    out_dir="03_fastqc2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q \
      -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc2:
  input:
    expand("03_fastqc2/{sample}_{rep}_trimmed_fastqc.html",sample=SAMPLE,rep=REP),
    expand("03_fastqc2/{sample}_{rep}_trimmed_fastqc.zip",sample=SAMPLE,rep=REP)
  output:
    "03_fastqc2/multiqc_report.html"
  log:
    "logs/fastqc2/multiqc2.log"
  threads: 1
  params:
    out_dir="03_fastqc2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o {params.out_dir} > {log} 2>&1"
    



