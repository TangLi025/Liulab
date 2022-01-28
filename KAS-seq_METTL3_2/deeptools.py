TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
#SAMPLE=["DMSO"]
SAMPLE=["DMSO","DRB","TRIP"]
MERGE=["","_merge"]

GENOME_2bit="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.2bit"
GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.gtf"
GTF_protein_coding="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.gene_coding.gtf"
BLACKLIST="/disk1/home/user_09/reference/annotation/hg19/hg19_blacklist.bed"

rule all:
  input:
    expand("07_deeptools/bamCoverage/{sample}_{treatment}_{rep}.bw",sample=SAMPLE,treatment="IP",rep=["merge"]),
    #expand("07_deeptools/plotFingerprint/{sample}_KAS-seq_plotFingerprint_{dup}.png",sample=SAMPLE,dup=DUP),
    #expand("07_deeptools/computeGCBias/{sample}_{treatment}_{rep}_{dup}.txt",sample=SAMPLE,treatment=TREATMENT,rep=REP,dup=DUP),
    #expand("07_deeptools/multiBigwigSummary/{sample}.npz",sample=SAMPLE),
    #expand("07_deeptools/plotCorrelation/{sample}.png",sample=SAMPLE),
    #expand("07_deeptools/plotCorrelation/{sample}.tab",sample=SAMPLE),
    #expand("07_deeptools/computeMatrix/{sample}.mat.gz",sample=SAMPLE),
    #expand("07_deeptools/computeMatrix/{sample}.tab",sample=SAMPLE),
    #expand("07_deeptools/computeMatrix/{sample}.bed",sample=SAMPLE),
    #expand("07_deeptools/plotProfile/{sample}.png",sample=SAMPLE),
    #expand("07_deeptools/plotHeatmap/{sample}.png",sample=SAMPLE),
    #"07_deeptools/plotHeatmap/compare.png",
    #"07_deeptools/plotProfile/compare_input.png"
    
    
rule faToTwoBit:
  input:
    GENOME
  output:
    GENOME_2bit
  shell:
    "/disk1/home/user_09/anaconda3/bin/faToTwoBit {input} {output}"

rule bam_merge:
  input:
    "04_bam_rmdup/{sample}_{treatment}_rep1.bam",
    "04_bam_rmdup/{sample}_{treatment}_rep2.bam"
  output:
    "04_bam_rmdup/{sample}_{treatment}_merge.bam"
  log:
    "logs/bam_merge/{sample}_{treatment}.log"
  threads:2
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools merge -f -@ {threads} {output} {input[0]} {input[1]} > {log} 2>&1"

rule bam_index:
  input:
    "04_bam_rmdup/{sample}_{treatment}_merge.bam"
  output:
    "04_bam_rmdup/{sample}_{treatment}_merge.bam.bai"
  log:
    "logs/bam_index_merge/{sample}_{treatment}_merge.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
 
rule computeGCBias:
  input:
    bam="04_bam_{dup}/{sample}_{treatment}_{rep}.bam",
    genome=GENOME_2bit
  output:
    freq="07_deeptools/computeGCBias/{sample}_{treatment}_{rep}_{dup}.txt",
    png="07_deeptools/computeGCBias/{sample}_{treatment}_{rep}_{dup}.png"
  log:
    "logs/computeGCBias/{sample}_{treatment}_{rep}_{dup}.log"
  params:
    genome=GENOME_2bit
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeGCBias \
      -b {input.bam} \
      --effectiveGenomeSize 2862010578 \
      --genome {input.genome} \
      --GCbiasFrequenciesFile {output.freq} \
      --fragmentLength 150 \
      --biasPlot {output.png} --plotFileFormat png \
      --regionSize 300 \
      --numberOfProcessors {threads} > {log} 2>&1"
      
rule multiBigwigSummary:
  input:
    "07_deeptools/bamCoverage/{sample}_IP_rep1.bw",
    "07_deeptools/bamCoverage/{sample}_IP_rep2.bw"
  output:
    npz="07_deeptools/multiBigwigSummary/{sample}.npz"
  params:
    blacklist=BLACKLIST
  log:
    "logs/multiBigwigSummary/{sample}.log"
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/multiBigwigSummary bins \
      -b {input} \
      --outFileName {output.npz} \
      --binSize 10000 \
      --numberOfProcessors {threads} \
      --blackListFileName {params.blacklist} \
       > {log} 2>&1"

rule plotCorrelation:
  input:
    npz="07_deeptools/multiBigwigSummary/{sample}.npz"
  output:
    png="07_deeptools/plotCorrelation/{sample}.png",
    tab="07_deeptools/plotCorrelation/{sample}.tab"
  log:
    "logs/plotCorrelation/{sample}.log"
  params:
    labels=["rep1","rep2"]
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotCorrelation \
      --corData {input.npz} \
      --corMethod pearson \
      --whatToPlot scatterplot \
      --plotFile {output.png} \
      --labels {params.labels} \
      --outFileCorMatrix {output.tab} \
      --xRange 0 75 \
      --yRange 0 75 \
      --removeOutliers \
       > {log} 2>&1"
       
rule computeMatrix_distribution:
  input:
    bw=expand("07_deeptools/bamCoverage/{sample}_{treatment}_{rep}.bw",sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="07_deeptools/computeMatrix/{sample}.mat.gz",
    tab="07_deeptools/computeMatrix/{sample}.tab",
    bed="07_deeptools/computeMatrix/{sample}.bed"
  log:
    "logs/computeMatrix_distribution/{sample}.log"
  params:
    gtf=GTF,
    blacklist=BLACKLIST,
    labels=["input_rep1","input_rep2","IP_rep1","IP_rep2"]
  threads: 45
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
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
      --samplesLabel {params.labels} \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile:
  input:
    mat="07_deeptools/computeMatrix/{sample}.mat.gz"
  output:
    png="07_deeptools/plotProfile/{sample}.png"
  log:
    "logs/plotProfile/{sample}.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"

rule plotHeatmap:
  input:
    mat="07_deeptools/computeMatrix/{sample}.mat.gz"
  output:
    png="07_deeptools/plotHeatmap/{sample}.png"
  log:
    "logs/plotHeatmap/{sample}.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotHeatmap \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"

rule computeMatrix_distribution_compare:
  input:
    bw=expand("07_deeptools/bamCoverage/{sample}_IP_merge.bw",sample=SAMPLE)
  output:
    mat="07_deeptools/computeMatrix/compare.mat.gz",
    tab="07_deeptools/computeMatrix/compare.tab",
    bed="07_deeptools/computeMatrix/compare.bed"
  log:
    "logs/computeMatrix_distribution/compare.log"
  params:
    gtf=GTF,
    blacklist=BLACKLIST,
    labels=["Native","DRB","TRIP"]
  threads: 45
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
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
      --samplesLabel {params.labels} \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotHeatmap_compare:
  input:
    mat="07_deeptools/computeMatrix/compare.mat.gz"
  output:
    png="07_deeptools/plotHeatmap/compare.png"
  log:
    "logs/plotHeatmap/compare.log"
  params:
    genes="genes",
    what_to_show="'heatmap and colorbar'"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotHeatmap \
      --sortUsingSamples 1 \
      --whatToShow {params.what_to_show} \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
      --colorMap 'Blues' \
       > {log} 2>&1"
       
rule computeMatrix_distribution_compare_input:
  input:
    expand("07_deeptools/bamCoverage/{sample}_IP_merge.bw",sample=SAMPLE),
    "07_deeptools/bamCoverage/DMSO_input_merge.bw"
  output:
    mat="07_deeptools/computeMatrix/compare_input.mat.gz",
    tab="07_deeptools/computeMatrix/compare_input.tab",
    bed="07_deeptools/computeMatrix/compare_input.bed"
  log:
    "logs/computeMatrix_distribution/compare_input.log"
  params:
    gtf=GTF_protein_coding,
    blacklist=BLACKLIST,
    labels=["Native","DRB","TRIP","input"]
  threads: 45
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
      --scoreFileName {input} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --regionBodyLength 6000 \
      --beforeRegionStartLength 3000  \
      --afterRegionStartLength 3000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --samplesLabel {params.labels} \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
      --metagene \
       > {log} 2>&1"      
       
rule plotProfile_compare_input:
  input:
    mat="07_deeptools/computeMatrix/compare_input.mat.gz"
  output:
    png="07_deeptools/plotProfile/compare_input.png"
  log:
    "logs/plotProfile/compare_input.log"
  params:
    genes="genes"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"
