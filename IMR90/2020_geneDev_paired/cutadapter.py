SAMPLE=["CTRL","SENE"]

REP=["rep1","rep2","rep3"]
READ=["1","2"]

DUP=["raw"]

INDEX="/disk1/home/user_09/reference/index/hisat2/hg38/hg38"

rule all:
  input:
    #expand("02_trim_galore/{sample}_{rep}_1_val_1.fq.gz",sample=SAMPLE,rep=REP),
    #expand("02_trim_galore/{sample}_{rep}_2_val_2.fq.gz",sample=SAMPLE,rep=REP),
    expand("04_bam_{dup}/{sample}_{rep}.bam.bai",dup=DUP,sample=SAMPLE,rep=REP),
    "01_fastqc1/multiqc_report.html",
    "03_fastqc2/multiqc_report.html"
    
rule fastqc1:
  input:
    "00_raw_fastq/{sample}_{rep}_{read}.fastq"
  output:
    "01_fastqc1/{sample}_{rep}_{read}_fastqc.html",
    "01_fastqc1/{sample}_{rep}_{read}_fastqc.zip"
  log:
    "logs/fastqc1/{sample}_{rep}_{read}.log"
  threads: 1
  params:
    out_dir="01_fastqc1"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc1:
  input:
    expand("01_fastqc1/{sample}_{rep}_{read}_fastqc.html", sample=SAMPLE,rep=REP,read=READ),
    expand("01_fastqc1/{sample}_{rep}_{read}_fastqc.zip", sample=SAMPLE,rep=REP,read=READ)
  output:
    "01_fastqc1/multiqc_report.html"
  log:
    "logs/fastqc1/multiqc1.log"
  threads: 1
  params:
    out_dir="01_fastqc1"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} -o {params.out_dir} > {log} 2>&1"
    
rule trim_galore:
  input:
    "00_raw_fastq/{sample}_{rep}_1.fastq",
    "00_raw_fastq/{sample}_{rep}_2.fastq"
  output:
    temp("02_trim_galore/{sample}_{rep}_1_val_1.fq"),
    temp("02_trim_galore/{sample}_{rep}_2_val_2.fq")
  params:
    output_dir="02_trim_galore"
  log:
    "logs/trim_galore/{sample}_{rep}.log"
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/trim-galore/bin/trim_galore \
      -o {params.output_dir} \
      --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt \
      -j {threads} --paired --trim1 {input} > {log} 2>&1"
      
rule fastqc2:
  input:
    "02_trim_galore/{sample}_{rep}_{read}_val_{read}.fq"
  output:
    "03_fastqc2/{sample}_{rep}_{read}_val_{read}_fastqc.html",
    "03_fastqc2/{sample}_{rep}_{read}_val_{read}_fastqc.zip"
  log:
    "logs/fastqc2/{sample}_{rep}_{read}.log"
  threads: 1
  params:
    out_dir="03_fastqc2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc2:
  input:
    expand("03_fastqc2/{sample}_{rep}_{read}_val_{read}_fastqc.html",sample=SAMPLE,rep=REP,read=READ),
    expand("03_fastqc2/{sample}_{rep}_{read}_val_{read}_fastqc.zip",sample=SAMPLE,rep=REP,read=READ)
  output:
    "03_fastqc2/multiqc_report.html"
  log:
    "logs/fastqc2/multiqc2.log"
  threads: 1
  params:
    out_dir="03_fastqc2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} -o {params.out_dir} > {log} 2>&1"

rule hisat2_mapping:
  input:
    "02_trim_galore/{sample}_{rep}_1_val_1.fq",
    "02_trim_galore/{sample}_{rep}_2_val_2.fq"
  output:
    bam="04_bam_raw/{sample}_{rep}.bam",
    summary="04_bam_raw/{sample}_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/hisat2_mapping/{sample}_{rep}.log"
  threads: 20
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/hisat2 -x {params.index} --rna-strandness RF \
    --summary-file {output.summary} \
    -p {threads} -1 {input[0]} -2 {input[1]} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -bS \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output.bam} 
    """
    
rule bam_index:
  input:
    "04_bam_{dup}/{sample}_{rep}.bam"
  output:
    "04_bam_{dup}/{sample}_{rep}.bam.bai"
  log:
    "logs/bam_index/{sample}_{rep}_{dup}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    
