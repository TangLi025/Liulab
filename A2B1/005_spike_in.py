SAMPLE=["Lysate","Result"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

rule all:
  input:
    # bowtie2比对并建立索引
    expand("12_spike_in/{sample}_{treatment}_{rep}.bam.bai",sample=SAMPLE,treatment=TREATMENT,rep=REP)
    
rule bowtie2_mapping_spike_in:
  input:
    "02_trim_galore/{sample}_{treatment}_{rep}_1_val_1.fq.gz",
    "02_trim_galore/{sample}_{treatment}_{rep}_2_val_2.fq.gz"
  output:
    temp("12_spike_in/{sample}_{treatment}_{rep}.sam"),
    "12_spike_in/{sample}_{treatment}_{rep}.bam",
    "12_spike_in/{sample}_{treatment}_{rep}_summary.txt"
  threads: 10
  log:
    "logs/12_spike_in/{sample}_{treatment}_{rep}.log"
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/bowtie2/bin/bowtie2 --no-unal -p {threads} \
    -x /disk1/home/user_09/reference/index/bowtie2/spike_in/spike_in -1 {input[0]} -2 {input[1]} \
     -S {output[0]} 1> {output[2]} 2>{log}
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -bS {output[0]}\
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output[1]} >> {log} 2>&1
    """
  
rule bam_index:
  input:
    "12_spike_in/{sample}_{treatment}_{rep}.bam"
  output:
    "12_spike_in/{sample}_{treatment}_{rep}.bam.bai"
  log:
    "logs/spike_in/bam_index_{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
  
