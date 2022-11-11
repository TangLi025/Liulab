SAMPLE=["CTRL","M3KO"]
TREATMENT=["input"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw"]
STRAND=["pos","neg"]

INDEX="/disk1/home/user_09/reference/index/hisat2/mm19_hisat2_ss_exon/mm19"

rule all:
  input:
    #expand("02_trim_galore/{sample}_{treatment}_{rep}_1_val_1.fq.gz",sample=SAMPLE,treatment=TREATMENT,rep=REP),
    #expand("02_trim_galore/{sample}_{treatment}_{rep}_2_val_2.fq.gz",sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("04_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai",dup=DUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
    #"01_fastqc1/multiqc_report.html",
    #"03_fastqc2/multiqc_report.html"
    
rule fastqc1:
  input:
    "00_raw_fastq/{sample}_{treatment}_{rep}_{read}.fq.gz"
  output:
    "01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.html",
    "01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.zip"
  log:
    "logs/fastqc1/{sample}_{treatment}_{rep}_{read}.log"
  threads: 1
  params:
    out_dir="01_fastqc1"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc1:
  input:
    expand("01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.html", sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ),
    expand("01_fastqc1/{sample}_{treatment}_{rep}_{read}_fastqc.zip", sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ)
  output:
    "01_fastqc1/multiqc_report.html"
  log:
    "logs/01_fastqc1/multiqc1.log"
  threads: 1
  params:
    out_dir="01_fastqc1"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/multiqc {input} -o {params.out_dir} > {log} 2>&1"
    
rule trim_galore:
  input:
    "00_raw_fastq/{sample}_{treatment}_{rep}_1.fastq",
    "00_raw_fastq/{sample}_{treatment}_{rep}_2.fastq"
  output:
    temp("02_trim_galore/{sample}_{treatment}_{rep}_1_val_1.fq"),
    temp("02_trim_galore/{sample}_{treatment}_{rep}_2_val_2.fq")
  params:
    output_dir="02_trim_galore"
  log:
    "logs/trim_galore/{sample}_{treatment}_{rep}.log"
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/trim-galore/bin/trim_galore \
      -o {params.output_dir} \
      --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt \
      -j {threads} --paired --clip_R1 10 --clip_R2 10 --trim1 {input} > {log} 2>&1"
      
rule fastqc2:
  input:
    "02_trim_galore/{sample}_{treatment}_{rep}_{read}_val_{read}.fq.gz"
  output:
    "03_fastqc2/{sample}_{treatment}_{rep}_{read}_val_{read}_fastqc.html",
    "03_fastqc2/{sample}_{treatment}_{rep}_{read}_val_{read}_fastqc.zip"
  log:
    "logs/fastqc2/{sample}_{treatment}_{rep}_{read}.log"
  threads: 1
  params:
    out_dir="03_fastqc2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/fastqc -t {threads} -q -o {params.out_dir} {input} > {log} 2>&1"
    
rule multiqc2:
  input:
    expand("03_fastqc2/{sample}_{treatment}_{rep}_{read}_val_{read}_fastqc.html",sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ),
    expand("03_fastqc2/{sample}_{treatment}_{rep}_{read}_val_{read}_fastqc.zip",sample=SAMPLE,treatment=TREATMENT,rep=REP,read=READ)
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
    "02_trim_galore/{sample}_{treatment}_{rep}_1_val_1.fq",
    "02_trim_galore/{sample}_{treatment}_{rep}_2_val_2.fq"
  output:
    bam="04_bam_raw/{sample}_{treatment}_{rep}.bam",
    summary="04_bam_raw/{sample}_{treatment}_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/hisat2_mapping/{sample}_{treatment}_{rep}.log"
  threads: 40
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/hisat2 -x {params.index} --rna-strandness RF \
    --summary-file {output.summary} \
    -p {threads} -1 {input[0]} -2 {input[1]} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -bS \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output.bam} 
    """
    
rule deduplicate:
  input:
    "04_bam_raw/{sample}_{treatment}_{rep}.bam"
  output:
    "04_bam_dedup/dedup_record/{sample}_{treatment}_{rep}.log",
    "04_bam_dedup/{sample}_{treatment}_{rep}.bam"
  log:
    "logs/bam_dedup/deduplicate_{sample}_{treatment}_{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/picard MarkDuplicates REMOVE_DUPLICATES=true \
      I={input} M={output[0]} O={output[1]} > {log} 2>&1"

rule bam_index:
  input:
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam"
  output:
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai"
  log:
    "logs/bam_index/{sample}_{treatment}_{rep}_{dup}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    
