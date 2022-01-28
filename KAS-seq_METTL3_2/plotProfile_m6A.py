GROUP=["METTL3_3","METTL3_2"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"

M6A_PLUS="/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/m6A_plus.bed"
M6A_MINUS="/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/m6A_minus.bed"
M6A_PEAK=expand("/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/Mettl3_{sample}_commen_peaks_mm19.narrowPeak",sample=["Control","KO1","KO2"])
  
BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

rule all:
  input:
    #expand("{group}/07_deeptools/computeMatrix_m6A/{group}.mat.gz",group=GROUP),
    #expand("{group}/07_deeptools/plotProfile_m6A/{group}.png",group=GROUP),
    #expand("{group}/07_deeptools/plotHeatmap_m6A/{group}.png",group=GROUP),
    #expand("{group}/07_deeptools/plotPCA_m6A/{group}.png",group=GROUP),
    #expand("{group}/07_deeptools/plotProfile_m6A/{group}_{type}.png",group=GROUP,type=["TSS","TES"]),
    expand("{group}/07_deeptools/plotProfile_m6A/m6A_peak_{group}.png",group=GROUP),
    
    
rule computeMatrix_distribution_m6A_peak:
  input:
    bw=expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix_m6A/m6A_peak_{group}.mat.gz",
    tab="{group}/07_deeptools/computeMatrix_m6A/m6A_peak_{group}.tab",
    bed="{group}/07_deeptools/computeMatrix_m6A/m6A_peak_{group}.bed"
  log:
    "{group}/logs/computeMatrix_distribution_m6A/m6A_peak_{group}.log"
  params:
    region=M6A_PEAK,
    blacklist=BLACKLIST,
  threads: 25
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix reference-point \
      --regionsFileName {params.region} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --referencePoint center \
      --beforeRegionStartLength 5000  \
      --afterRegionStartLength 5000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --smartLabels \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_m6A_peak:
  input:
    mat="{group}/07_deeptools/computeMatrix_m6A/m6A_peak_{group}.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile_m6A/m6A_peak_{group}.png"
  log:
    "{group}/logs/plotProfile_m6A/m6A_peak_{group}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"
    
rule computeMatrix_distribution_TSS_TES:
  input:
    bw=expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_IP_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix_m6A/{group}_{type}.mat.gz",
    tab="{group}/07_deeptools/computeMatrix_m6A/{group}_{type}.tab",
    bed="{group}/07_deeptools/computeMatrix_m6A/{group}_{type}.bed"
  log:
    "{group}/logs/computeMatrix_distribution_m6A/{group}_{type}.log"
  params:
    type=r'{type}',
    region=[M6A_PLUS,M6A_MINUS],
    blacklist=BLACKLIST,
  threads: 20
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix reference-point \
      --regionsFileName {params.region} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --referencePoint {params.type} \
      --beforeRegionStartLength 5000  \
      --afterRegionStartLength 5000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --smartLabels \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_TSS_TES:
  input:
    mat="{group}/07_deeptools/computeMatrix_m6A/{group}_{type}.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile_m6A/{group}_{type}.png"
  log:
    "{group}/logs/plotProfile_m6A/{group}_{type}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"   

rule computeMatrix_distribution:
  input:
    bw=expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix_m6A/{group}.mat.gz",
    tab="{group}/07_deeptools/computeMatrix_m6A/{group}.tab",
    bed="{group}/07_deeptools/computeMatrix_m6A/{group}.bed"
  log:
    "{group}/logs/computeMatrix_distribution_m6A/{group}.log"
  params:
    region=[M6A_PLUS,M6A_MINUS],
    blacklist=BLACKLIST,
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.region} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --regionBodyLength 6000 \
      --beforeRegionStartLength 3000  \
      --afterRegionStartLength 3000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --smartLabels \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile:
  input:
    mat="{group}/07_deeptools/computeMatrix_m6A/{group}.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile_m6A/{group}.png"
  log:
    "{group}/logs/plotProfile_m6A/{group}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"
       
rule plotHeatmap:
  input:
    mat="{group}/07_deeptools/computeMatrix_m6A/{group}.mat.gz"
  output:
    png="{group}/07_deeptools/plotHeatmap_m6A/{group}.png"
  log:
    "{group}/logs/plotHeatmap_m6A/{group}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotHeatmap \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
      --colorMap 'Blues' \
       > {log} 2>&1"

rule multiBigwigSummary:
  input:
    expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=r"{group}",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    npz="{group}/07_deeptools/multiBigwigSummary_m6A/{group}.npz",
    tab="{group}/07_deeptools/multiBigwigSummary_m6A/{group}.tab"
  params:
    blacklist=BLACKLIST
  log:
    "{group}/logs/multiBigwigSummary_m6A/{group}.log"
  threads: 5
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/multiBigwigSummary bins \
      -b {input} \
      --outFileName {output.npz} \
      --outRawCounts {output.tab} \
      --binSize 10000 \
      --numberOfProcessors {threads} \
      --blackListFileName {params.blacklist} \
       > {log} 2>&1"

rule plotPCA:
  input:
    "{group}/07_deeptools/multiBigwigSummary_m6A/{group}.npz"
  output:
    "{group}/07_deeptools/plotPCA_m6A/{group}.png"
  log:
    "{group}/logs/plotPCA_m6A/{group}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotPCA \
      -in {input} \
      --plotHeight 15 --plotWidth 10 \
      --plotFile {output} \
      > {log} 2>&1"
