GROUP=["KAS-seq_ALKBH3","KAS-seq_ALKBH5"]
SAMPLE=["Ctrl","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]


rule all:
  input:
    expand("{group}/02_trim_galore/{group}_{sample}_{treatment}_{rep}_trimmed.fq.gz",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)

'''
rule fasterq_dump:
  input:
    "{group}/00_raw_sra/{sample}_{treatment}_{rep}.sra"
  output:
    temp("{group}/00_raw_fastq/{sample}_{treatment}_{rep}.fastq")
  log:
    "{group}/logs/fasterq_dump/{sample}_{treatment}_{rep}.log"
  threads:10
  shell:
     "/disk1/home/user_09/anaconda3/envs/m6A/bin/fasterq-dump -e {threads} --split-3 -O {wildcards.group}/00_raw_fastq {input} > {log} 2>&1"
'''
 
rule fastqc1:
  input:
    "{group}/00_raw_fastq/{group}_{sample}_{treatment}_{rep}.fastq.gz"
  output:
    "{group}/01_fastqc1/{group}_{sample}_{treatment}_{rep}_fastqc.html",
    "{group}/01_fastqc1/{group}_{sample}_{treatment}_{rep}_fastqc.zip"
  log:
    "{group}/logs/fastqc1/{group}_{sample}_{treatment}_{rep}.log"
  params:
    out_dir="{group}/01_fastqc1"
  threads: 2
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc1:
  input:
    expand("{group}/01_fastqc1/{group}_{sample}_{treatment}_{rep}_fastqc.html",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("{group}/01_fastqc1/{group}_{sample}_{treatment}_{rep}_fastqc.zip",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    "{group}/01_fastqc1/multiqc_report.html"
  log:
    "{group}/logs/fastqc1/multiqc1.log"
  threads: 1
  params:
    out_dir="{group}/01_fastqc1"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o {params.out_dir} > {log} 2>&1"

rule trim_galore:
  input:
    "{group}/00_raw_fastq/{group}_{sample}_{treatment}_{rep}.fastq.gz"
  output:
    "{group}/02_trim_galore/{group}_{sample}_{treatment}_{rep}_trimmed.fq.gz"
  params:
    output_dir="{group}/02_trim_galore"
  log:
    "{group}/logs/trim_galore/{group}_{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/trim-galore/bin/trim_galore \
    --gzip \
    --length 50 -j {threads} \
    --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt \
    --output_dir {params.output_dir} {input} > {log} 2>&1"

rule fastqc2:
  input:
    "{group}/02_trim_galore/{group}_{sample}_{treatment}_{rep}_trimmed.fq.gz"
  output:
    "{group}/03_fastqc2/{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.html",
    "{group}/03_fastqc2/{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.zip"
  log:
    "{group}/logs/fastqc2/{group}_{sample}_{treatment}_{rep}.log"
  threads: 2
  params:
    out_dir="{group}/03_fastqc2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q \
      -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc2:
  input:
    expand("{group}/03_fastqc2/{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.html",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("{group}/03_fastqc2/{group}_{sample}_{treatment}_{rep}_trimmed_fastqc.zip",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    "{group}/03_fastqc2/multiqc_report.html"
  log:
    "{group}/logs/fastqc2/multiqc2.log"
  threads: 1
  params:
    out_dir="{group}/03_fastqc2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o {params.out_dir} > {log} 2>&1"
    



