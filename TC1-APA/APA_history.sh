

##APA trim-galore

for i in {1..10}; do trim_galore -o 02_trim_galore/ --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt -j 10 --paired 00_raw_fastq/TC1-AP${i}_L1_1.fq.gz 00_raw_fastq/TC1-AP${i}_L1_2.fq.gz ; done &







 for i in {1..10}; do trim_galore -a A{20} -a2 T{20} --path_to_cutadapt /disk1/home/user_09/anaconda3/envs/trim-galore/bin/cutadapt -j 15 --polyA --length 15 --paired 02_trim_galore/TC1-AP${i}_L1_1_val_1.fq.gz 02_trim_galore/TC1-AP${i}_L1_2_val_2.fq.gz; done &
 
 for i in {1..10}; do zcat TC1-AP${i}_L1_1_val_1_val_1.fq.gz | grep -A 3 PolyA | grep -v ^-- > TC1-AP${i}_1_polyA_trimmed.fastq; done &
 
 for i in {1..10}; do zcat TC1-AP${i}_L1_2_val_2_val_2.fq.gz | grep -A 3 PolyA | grep -v ^-- > TC1-AP${i}_2_polyA_trimmed.fastq; done &
 
 for i in {1..10}; do fastp -w 10 -i TC1-AP${i}_1_polyA_trimmed.fastq -o TC1-AP${i}_1_polyA_trimmed_paired.fastq -I TC1-AP${i}_2_polyA_trimmed.fastq -O TC1-AP${i}_2_polyA_trimmed_paired.fastq; done &