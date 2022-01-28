GROUP=["METTL3"]
SAMPLE=["CTRL","KO1","KO2"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

rule all:
  input:
    expand("00_raw_fastq/{group}_{sample}_{treatment}_{rep}_1.fastq",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep-REP)
    #"01_fastqc1/multiqc1/multiqc_report.html",
    #expand("02_trim_galore/{sample}_{treatment}_{rep}_{read}.fastq",sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ),
    #"03_fastqc2/multiqc2/multiqc_report.html"

rule fasterq_dump:
  input:
    "00_raw_sra/{group}_{sample}_{treatment}_{rep}.sra"
  output:
    "00_raw_fastq/{group}_{sample}_{treatment}_{rep}_1.fastq",
    "00_raw_fastq/{group}_{sample}_{treatment}_{rep}_2.fastq"
  log:
    "logs/fasterq_dump/{sample}_{treatment}_{rep}.log"
  threads:10
  shell:
     "/disk1/home/user_09/anaconda3/envs/m6A/bin/fasterq-dump -e {threads} --split-3 -O 00_raw_fastq {input} > {log} 2>&1"

"""  
rule fastqc1:
  input:
    "00_raw_fastq/{sample}_{treatment}_{rep}_{read}.fq.gz"
  output:
    "01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.html",
    "01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.zip"
  log:
    "logs/fastqc1/{sample}_{treatment}_{rep}_{read}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o 01_fastqc1 {input} > {log} 2>&1"

rule multiqc1:
  input:
    expand("01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.html",sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ),
    expand("01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.zip",sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ)
  output:
    "01_fastqc1/multiqc1/multiqc_report.html"
  log:
    "logs/fastqc1/multiqc1.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o 01_fastqc1/multiqc1/ > {log} 2>&1"

rule fastp:
  input:
    "00_raw_fastq/{sample}_{treatment}_{rep}_1.fq.gz",
    "00_raw_fastq/{sample}_{treatment}_{rep}_2.fq.gz"
  output:
    "02_fastp/{sample}_{treatment}_{rep}_1.fastq",
    "02_fastp/{sample}_{treatment}_{rep}_2.fastq"
  log:
    "logs/02_fastp/{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastp -w {threads} -h ./02_fastp/{wildcards.sample}.html -i {input[0]} -o {output[0]} -I {input[1]} -O {output[1]} > {log} 2>&1"

rule fastqc2:
  input:
    "02_trim_galore/{sample}_{treatment}_{rep}_{read}_trimmed.fq"
  output:
    "03_fastqc2/{sample}_{treatment}_{rep}_{read}_trimmed_fastqc.html",
    "03_fastqc2/{sample}_{treatment}_{rep}_{read}_trimmed_fastqc.zip"
  log:
    "logs/fastqc2/{sample}_{treatment}_{rep}_{read}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q \
      -o 03_fastqc2 {input} > {log} 2>&1"
    
rule multiqc2:
  input:
    expand("03_fastqc2/{sample}_{treatment}_{rep}_{read}_trimmed_fastqc.html",sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ),
    expand("03_fastqc2/{sample}_{treatment}_{rep}_{read}_trimmed_fastqc.zip",sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ)
  output:
    "03_fastqc2/multiqc2/multiqc_report.html"
  log:
    "logs/fastqc2/multiqc2.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} \
      -o 03_fastqc2/multiqc2 > {log} 2>&1"
"""
