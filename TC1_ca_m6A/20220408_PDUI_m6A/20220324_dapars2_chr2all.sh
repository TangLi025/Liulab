#!/bin/bash

for i in {1..22}
do
  sed -i '1d' /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_chr${i}/DaPars2_mRNA_result_temp.chr${i}.txt
  cat /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_chr${i}/DaPars2_mRNA_result_temp.chr${i}.txt >> \
    /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_chrX/DaPars2_mRNA_result_temp.chrX.txt
done

sed -i '1d' /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_chrY/DaPars2_mRNA_result_temp.chrY.txt

cat /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_chrY/DaPars2_mRNA_result_temp.chrY.txt >> \
  /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_chrX/DaPars2_mRNA_result_temp.chrX.txt
  
cp /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_chrX/DaPars2_mRNA_result_temp.chrX.txt \
  /disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_mRNA_new/DaPars2_mRNA_result.txt