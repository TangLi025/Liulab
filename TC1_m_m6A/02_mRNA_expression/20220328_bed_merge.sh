#!/bin/sh
cd ~/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/
#cd ~/Data/TC1-planB/05_bam_change_index_hisat2/02_bam_dedup/
mkdir ./10_bed_merge
## part 1: replicates intersect (for peak merge between samples and featurecount)
mkdir ./10_bed_merge/00_common_peaks
for i in {"p0","p5","p10","rp2"}
do
echo 08_bed_filtered/${i}_rep1_peaks.bed
echo 08_bed_filtered/${i}_rep2_peaks.bed
bedtools intersect -a 08_bed_filtered/${i}_rep1_peaks.bed -b 08_bed_filtered/${i}_rep2_peaks.bed -f 0.5 -F 0.5 -e -s > 10_bed_merge/00_common_peaks/${i}_rep1_rep2_common_peaks.bed
echo 10_bed_merge/00_common_peaks/${i}_rep1_rep2_common_peaks.bed
cat 10_bed_merge/00_common_peaks/${i}_rep1_rep2_common_peaks.bed | wc -l
done


## part 2: merge between replicates (for reproducibility)

mkdir ./10_bed_merge/01_union_peaks
for i in {"p0","p5","p10","rp2"}
do
cat ./08_bed_filtered/${i}_rep1_peaks.bed ./08_bed_filtered/${i}_rep2_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct | awk -F "\t" -v OFS="\t" '{print "peak"NR,$1,$2,$3,$4}'> ./10_bed_merge/01_union_peaks/${i}.union.saf
done

## part 3: merge between samples (for diff peak)

cd ./10_bed_merge/00_common_peaks
mkdir ../02_merge_peaks
cat p0_rep1_rep2_common_peaks.bed p5_rep1_rep2_common_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct | awk -F "\t" -v OFS="\t" '{print "peak"NR,$1,$2,$3,$4}'> ../02_merge_peaks/p0p5.merge.saf
cat p0_rep1_rep2_common_peaks.bed p10_rep1_rep2_common_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct | awk -F "\t" -v OFS="\t" '{print "peak"NR,$1,$2,$3,$4}'> ../02_merge_peaks/p0p10.merge.saf
cat rp2_rep1_rep2_common_peaks.bed p5_rep1_rep2_common_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct | awk -F "\t" -v OFS="\t" '{print "peak"NR,$1,$2,$3,$4}'> ../02_merge_peaks/rp2p5.merge.saf
cat rp2_rep1_rep2_common_peaks.bed p10_rep1_rep2_common_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct | awk -F "\t" -v OFS="\t" '{print "peak"NR,$1,$2,$3,$4}'> ../02_merge_peaks/rp2p10.merge.saf


