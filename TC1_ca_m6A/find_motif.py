SAMPLE=["TC1_P0","TC1_P5","TC1_P10","TC1_rP2"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","dedup"]
STRAND=["pos","neg"]

GENOME="/disk1/home/user_09/reference/genome/mm/GRCm39.genome.fa"

rule all:
  input:
    # homer call motif
    expand("09_motif/{dup}/{sample}_{rep}/homerResults.html",sample=SAMPLE,rep=REP,dup=DUP),
    
    # bed_merge
    expand("10_bed_merge/{dup}/{sample}_{rep}_unique_peaks.bed",rep=REP,sample=SAMPLE,dup=DUP),
    expand("10_bed_merge/{dup}/{sample}_{rep}_common_peaks.bed",rep=REP,sample=SAMPLE,dup=DUP),
    expand("10_bed_merge/{dup}/{sample}_common_peaks.bed",sample=SAMPLE,dup=DUP),
    
    # bam_merge
    expand("11_bam_merge/{dup}/{sample}_{treatment}_{rep}.bam",treatment=TREATMENT,rep=REP,sample=SAMPLE,dup=DUP),
    
rule find_motif:
  input:
    "08_bed_filtered/{dup}/{sample}_{rep}_peaks.bed"
  output:
    "09_motif/{dup}/{sample}_{rep}/homerResults.html"
  log:
    "logs/09_motif/{dup}/{sample}_{rep}_find_motif.log"
  threads:2
  params:
    genome=GENOME,
    out_dir="09_motif/{dup}/{sample}_{rep}"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/findMotifsGenome.pl {input} \
      {params.genome} {params.out_dir} \
      -rna -p {threads} -len 5,6,7 > {log} 2>&1"

rule bed_merge:
  input:
    "08_bed_filtered/{dup}/{sample}_rep1_peaks.bed",
    "08_bed_filtered/{dup}/{sample}_rep2_peaks.bed"
  output:
    "10_bed_merge/{dup}/{sample}_rep1_unique_peaks.bed",
    "10_bed_merge/{dup}/{sample}_rep2_unique_peaks.bed",
    "10_bed_merge/{dup}/{sample}_rep1_common_peaks.bed",
    "10_bed_merge/{dup}/{sample}_rep2_common_peaks.bed",
    "10_bed_merge/{dup}/{sample}_common_peaks.bed"
  log:
    "logs/10_bed_merge/{dup}/{sample}.log"
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[0]} -b {input[1]} -v -s \
    > {output[0]}
    echo "{output[0]}:" > {log} 2>&1
    cat {output[0]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[1]} -b {input[0]} -v -s \
    > {output[1]}
    echo "{output[1]}:" >> {log} 2>&1
    cat {output[1]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[0]} -b {output[0]} -v -s \
    > {output[2]}
    echo "{output[2]}:" >> {log} 2>&1
    cat {output[2]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[1]} -b {output[1]} -v -s \
    > {output[3]}
    echo "{output[3]}:" >> {log} 2>&1
    cat {output[3]} | wc -l >> {log} 2>&1
  
    cat {output[2]} {output[3]} \
    | sort -k1,1 -k2,2n \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools merge -s -c 6 -o distinct \
    | awk -F "\t" -v OFS="\t" '{{print $1,$2,$3,NR,".",$4}}' \
    > {output[4]}
    echo "{output[4]}:" >> {log} 2>&1
    cat {output[4]} | wc -l >> {log} 2>&1
    """

rule bam_merge:
  input:
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_pos.bam",
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_neg.bam"
  output:
    "11_bam_merge/{dup}/{sample}_{treatment}_{rep}.bam"
  log:
    "logs/11_bam_merge/{dup}/{sample}_{treatment}_{rep}.log"
  threads:4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools merge -f -@ {threads} {output} {input[0]} {input[1]} > {log} 2>&1"
