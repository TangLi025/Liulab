SAMPLE=["H3K27ac","H3K4me1"]
TREATMENT=["input","IP"]
REP=["rep"]

rule all:
  input:
    expand("06_macs2/{sample}_peaks.narrowPeak",sample=SAMPLE)
    
rule macs2_callpeak:
  input:
    "04_bam_unique/{sample}_IP_rep.bam",
    "04_bam_unique/CTRL_input_rep.bam"
  output:
    "06_macs2/{sample}_peaks.xls",
    "06_macs2/{sample}_peaks.narrowPeak"
  log:
    "logs/macs2_callpeak/{sample}.log"
  params:
    out_name="{sample}",
    out_dir="06_macs2"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n {params.out_name} \
      -f BAM --verbose 3 --nomodel --extsize 50 \
      -g mm -q 0.01 \
      --outdir {params.out_dir} > {log} 2>&1"
