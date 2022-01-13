GROUP=["METTL3_2"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

rule all:
  input:
    "07_deeptools/plotFingerprint/KAS-seq_plotFingerprint_rmdup.png"


rule plotFingerprint:
  input:
    expand("04_bam_rmdup/KAS-seq_{group}_{sample}_{treatment}_{rep}.bam",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    png="07_deeptools/plotFingerprint/KAS-seq_plotFingerprint_rmdup.png",
    tab="07_deeptools/plotFingerprint/KAS-seq_plotFingerprint_rmdup.tab"
  log:
    "logs/plotFingerprint/KAS-seq_plotFingerprint_rmdup.log"
  params:
    labels=expand("{samples}_{treatment}_{rep}",samples=SAMPLE,treatment=["input","KAS-seq"],rep=REP),
    title="""Fingerprints of KAS-seq data""",
  threads: 20
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/deeptools/bin/plotFingerprint \
      -b {input} --labels {params.labels} --minMappingQuality 30 \
      --skipZeros --region 1 --numberOfSamples 500000 \
      -T "KAS-seq Fingerprint plot" \
      --plotFile {output.png} --plotFileFormat png \
      --outRawCounts {output.tab} \
      --numberOfProcessors {threads} \
      > {log} 2>&1
    """
