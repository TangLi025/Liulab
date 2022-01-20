SAMPLE=["Lysate","Result"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","rmdup"]
STRAND=["pos","neg"]

rule all:
  input:
    expand("08_bed_filtered/{dup}/{sample}_{rep}_peaks.bed",dup=DUP,sample=SAMPLE,rep=REP)
    
rule bam_separate:
  input:
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam",
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai"
  output:
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_83.bam"),
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_163.bam"),
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_99.bam"),
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_147.bam"),
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_pos.bam",
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_neg.bam"
  threads: 4
  shell:
    """
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 83 {input[0]} 1> {output[0]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 163 {input[0]} 1> {output[1]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 99 {input[0]} 1> {output[2]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 147 {input[0]} 1> {output[3]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools merge -@ {threads} {output[4]} {output[0]} {output[1]} 
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools merge -@ {threads} {output[5]} {output[2]} {output[3]}
  	"""
    
rule bam_separated_index:
  input:
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_{strand}.bam"
  output:
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_{strand}.bam.bai"
  log:
    "logs/05_bam_{dup}_separated/bam_{strand}_index_{sample}_{treatment}_{rep}.log"
  threads: 2
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"

rule macs2_callpeak:
  input:
    "05_bam_{dup}_separated/{sample}_IP_{rep}_{strand}.bam",
    "05_bam_{dup}_separated/{sample}_input_{rep}_{strand}.bam"
  output:
    "06_macs2/{dup}/{sample}_{rep}_{strand}_summits.bed",
    "06_macs2/{dup}/{sample}_{rep}_{strand}_peaks.xls",
    "06_macs2/{dup}/{sample}_{rep}_{strand}_peaks.narrowPeak",
    "06_macs2/{dup}/{sample}_{rep}_{strand}_control_lambda.bdg",
    "06_macs2/{dup}/{sample}_{rep}_{strand}_treat_pileup.bdg"
  log:
    "logs/06_macs2/{dup}/{sample}_{rep}_{strand}.log"
  params:
    out_name="{sample}_{rep}_{strand}",
    out_dir="06_macs2/{dup}"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n {params.out_name} \
      -f BAM --verbose 3 --nomodel --extsize 150 \
      -g mm -q 0.01 \
      --outdir {params.out_dir} > {log} 2>&1"

rule bed_modify:
  input:
    "06_macs2/{dup}/{sample}_{rep}_neg_peaks.narrowPeak",
    "06_macs2/{dup}/{sample}_{rep}_pos_peaks.narrowPeak"
  output:
    "06_macs2/{dup}/bed_modify/{sample}_{rep}_neg_peaks.narrowPeak",
    "06_macs2/{dup}/bed_modify/{sample}_{rep}_pos_peaks.narrowPeak",
    "06_macs2/{dup}/bed_modify/{sample}_{rep}_peaks.narrowPeak",
    "06_macs2/{dup}/bed_modify/{sample}_{rep}_peaks.bed"
  log:
    "logs/06_macs2/{dup}/bed_modify/{sample}_{rep}.log"
  shell:
    """
    echo "primary peak number[neg]:" > {log} 2>&1
    cat {input[0]} | wc -l >> {log} 2>&1
    cat {input[0]} \
    | awk '$1 ~ /^chr[0-9]*$/ || $1 ~ /^chr[X|Y]$/' | awk '$7>=1' \
    | awk '$6="-"' | awk -v OFS="\t" '{{print $0}}' \
    1> {output[0]} 2>> {log}
    echo "modified peak number[neg]:" >> {log} 2>&1
    cat {output[0]} | wc -l >> {log} 2>&1
    
    echo "primary peak number[pos]:" >> {log} 2>&1
    cat {input[1]} | wc -l >> {log} 2>&1
    cat {input[1]} \
    | awk '$1 ~ /^chr[0-9]*$/ || $1 ~ /^chr[X|Y]$/' | awk '$7>=1' \
    | awk '$6="+"' | awk -v OFS="\t" '{{print $0}}' \
    1> {output[1]} 2>> {log}
    echo "modified peak number[pos]:" >> {log} 2>&1
    cat {output[1]} | wc -l >> {log} 2>&1
    
    cat {output[0]} {output[1]} \
    | sort -k1,1 -k2,2n | awk -v OFS="\t" '{{print $0}}' 1> {output[2]}
    cat {output[2]} | awk -v OFS="\t" '{{print $1,$2,$3,$4,$5,$6}}' > {output[3]}
    cat {output[3]} | wc -l >> {log} 2>&1
    """
    
rule bed_filter:
  input:
    "06_macs2/{dup}/bed_modify/{sample}_{rep}_peaks.bed",
    "05_bam_{dup}_separated/{sample}_input_{rep}_neg.bam",
    "05_bam_{dup}_separated/{sample}_input_{rep}_pos.bam"
  output:
    "07_featurecount_input_read_in_peak/{dup}/{sample}_{rep}.csv",
    "08_bed_filtered/{dup}/{sample}_{rep}_peaks.bed"
  log:
    "logs/08_bed_filtered/{dup}/{sample}_{rep}.log"
  script:
    "scripts/bed_filter.R"
    

      
