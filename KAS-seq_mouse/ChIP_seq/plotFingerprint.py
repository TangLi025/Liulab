SAMPLE=["CTRL_input","H3K27ac_IP","H3K4me1_IP"]
REP=["rep"]

rule all:
  input:
    "07_deeptools/plotFingerprint/plotFingerprint_unique.png"


rule plotFingerprint:
  input:
    expand("04_bam_unique/{sample}_{rep}.bam",sample=SAMPLE,rep=REP)
  output:
    png="07_deeptools/plotFingerprint/plotFingerprint_unique.png",
    tab="07_deeptools/plotFingerprint/plotFingerprint_unique.tab"
  log:
    "logs/plotFingerprint/plotFingerprint_unique.log"
  threads: 100
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/deeptools/bin/plotFingerprint \
      -b {input} --minMappingQuality 30 \
      --skipZeros --region 1 --numberOfSamples 500000 \
      -T "Fingerprints of ChIP-seq data" \
      --plotFile {output.png} --plotFileFormat png \
      --outRawCounts {output.tab} \
      --numberOfProcessors {threads} \
      > {log} 2>&1
    """
