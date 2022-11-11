#!/bin/bash

cd /disk/user_09/reference/annotation/mm19/repeats/repeats_family_separate

for family in `cat /disk/user_09/reference/annotation/mm19/repeats/family_name.txt`
do
    echo ${family}
    awk -v OFS='\t' '{print $11,$6,$7,$8,$10}' mm19_repeats_all_${family}.bed > mm19_repeats_${family}.saf
done