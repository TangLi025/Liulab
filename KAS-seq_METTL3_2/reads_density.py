GROUP=["METTL3_2","METTL3_3","ALKBH3","ALKBH5"]
SAMPLE=["CTRL","KO"]

TREATMENT=["IP","input"]
REP=["rep1","rep2"]


BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

BED_GB="/disk1/home/user_09/reference/annotation/mm19/tag_density_KAS/mm19_Refseq.GB.bed"
BED_TSS="/disk1/home/user_09/reference/annotation/mm19/tag_density_KAS/mm19_Refseq.TSS.bed"

rule all:
    input:
        expand("{group}/07_deeptools/reads_density/reads_density_TSS.tab",group=GROUP),
        expand("{group}/07_deeptools/reads_density/reads_density_GB.tab",group=GROUP)

rule reads_density_TSS:
    input:
        expand("{group}/05_bedtools/bigWig/{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
    output:
        npz="{group}/07_deeptools/reads_density/reads_density_TSS.npz",
        tab="{group}/07_deeptools/reads_density/reads_density_TSS.tab",
    params:
        bed=BED_TSS,
        bl=BLACKLIST
    log:
        "{group}/logs/reads_density/reads_density_TSS.log"
    threads: 50
    shell:
        "/disk1/home/user_09/anaconda3/envs/deeptools/bin/multiBigwigSummary \
            BED-file -b {input} -o {output.npz} \
            --BED {params.bed} --smartLabels \
            --blackListFileName {params.bl} \
            --numberOfProcessors {threads} \
            --outRawCounts {output.tab} > {log} 2>&1"
            
rule reads_density_GB:
    input:
        expand("{group}/05_bedtools/bigWig/{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
    output:
        npz="{group}/07_deeptools/reads_density/reads_density_GB.npz",
        tab="{group}/07_deeptools/reads_density/reads_density_GB.tab",
    params:
        bed=BED_GB,
        bl=BLACKLIST
    log:
        "{group}/logs/reads_density/reads_density_GB.log"
    threads: 50
    shell:
        "/disk1/home/user_09/anaconda3/envs/deeptools/bin/multiBigwigSummary \
            BED-file -b {input} -o {output.npz} \
            --BED {params.bed} --smartLabels \
            --blackListFileName {params.bl} \
            --numberOfProcessors {threads} \
            --outRawCounts {output.tab} > {log} 2>&1"
