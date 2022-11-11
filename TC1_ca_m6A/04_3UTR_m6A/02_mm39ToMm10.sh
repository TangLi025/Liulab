#!/bin/bash

for i in {p0p5,p0p10,rp2p5,rp2p10}
do
/disk/user_09/reference/annotation/liftOver/liftOver \
    /disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/02_merge_peaks/${i}.merge.bed \
    /disk/user_09/reference/annotation/liftOver/mm39ToMm10.over.chain.gz \
    /disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/02_merge_peaks/mm39ToMm10/${i}.merge.bed \
    /disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/02_merge_peaks/mm39ToMm10/${i}.merge.unmap.bed
done