GROUP=["mESC"]
SAMPLE=["Ctrl","DRB"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

rule all:
  input:
    expand("{group}/07_deeptools/plotFingerprint/KAS-seq_plotFingerprint_rmdup.png",group=GROUP)


rule plotFingerprint:
  input:
    expand("{group}/04_bam_unique/{group}_{sample}_{treatment}_{rep}.bam",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    png="{group}/07_deeptools/plotFingerprint/KAS-seq_plotFingerprint_rmdup.png",
    tab="{group}/07_deeptools/plotFingerprint/KAS-seq_plotFingerprint_rmdup.tab"
  log:
    "{group}/logs/plotFingerprint/KAS-seq_plotFingerprint_rmdup.log"
  params:
    labels=expand("{samples}_{treatment}_{rep}",samples=SAMPLE,treatment=["input","KAS-seq"],rep=REP),
    title=r"""Fingerprints of KAS-seq data""",
  threads: 20
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/deeptools/bin/plotFingerprint \
      -b {input} --labels {params.labels} --minMappingQuality 30 \
      --skipZeros --region 1 --numberOfSamples 500000 \
      -T "Fingerprints of KAS-seq data" \
      --plotFile {output.png} --plotFileFormat png \
      --outRawCounts {output.tab} \
      --numberOfProcessors {threads} \
      > {log} 2>&1
    """
