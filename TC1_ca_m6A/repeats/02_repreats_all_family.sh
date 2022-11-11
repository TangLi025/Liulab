#!/bin/bash

cd /disk/user_09/reference/annotation/mm19/repeats/
# mkdir repeats_family_separate
for family in `cat /disk/user_09/reference/annotation/mm19/repeats/family_name.txt`
do
    echo ${family}
    cat /disk/user_09/reference/annotation/mm19/repeats/mm19_repeats_all.bed |
        grep ${family}$'\t' > repeats_family_separate/mm19_repeats_all_${family}.bed
done