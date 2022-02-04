SAMPLE=["Lysate","Result"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]

DUP=["raw","dedup"]

rule all:
  input:
    expand("07_deeptools/plotFingerprint/plotFingerprint_{dup}.png",dup=DUP)

rule plotFingerprint:
  input:
    expand("04_bam_{dup}/{sample}_{treatment}_{rep}.bam",dup=r'{dup}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    png="07_deeptools/plotFingerprint/plotFingerprint_{dup}.png",
    tab="07_deeptools/plotFingerprint/plotFingerprint_{dup}.tab"
  log:
    "logs/plotFingerprint/plotFingerprint_{dup}.log"
  params:
    labels=expand("{samples}_{treatment}_{rep}",samples=SAMPLE,treatment=["input","m6A"],rep=REP),
    title=r"""Fingerprints of MeRIP data""",
  threads: 50
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/deeptools/bin/plotFingerprint \
      -b {input} --labels {params.labels} --minMappingQuality 30 \
      --skipZeros --region 1 --numberOfSamples 500000 \
      -T "Fingerprints of MeRIP data" \
      --plotFile {output.png} --plotFileFormat png \
      --outRawCounts {output.tab} \
      --numberOfProcessors {threads} \
      > {log} 2>&1
    """
