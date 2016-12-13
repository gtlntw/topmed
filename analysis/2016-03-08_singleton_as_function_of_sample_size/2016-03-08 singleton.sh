#!/bin/bash

#random  
split -l 1000 unrelated_degree_2.txt unrelated_degree_2_

#generate ac
PATH =$PATH\:/net/snowwhite/home/khlin/bin
cmd="sbatch --partition=topmed  --job-name=topmedsingleton -D $PWD -o $PWD/log/slurm_%j.out --wrap"
file_name=("aa" "ab" "ac" "ad" "ae" "af" "ag")
for x in {1..7}
do
	> unrelated_degree_2_${x}000.txt
	for y in `seq 1 $x`
	do
		idx=${y}-1
	 	cat unrelated_degree_2_${file_name[idx]} >> unrelated_degree_2_${x}000.txt
	done
done

for x in {1..7}
do
	for chr in {1..22}
	do
	 eval "$cmd \"bcftools view -f PASS -S unrelated_degree_2_${x}000.txt /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.genotypes.vcf.gz | bcftools query -i 'INFO/AC==1' -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' > temp/topmed_freeze2_${x}000.chr${chr}.overlap_removed.svm_pass.stats\"" 
	done
done

#count singleton
for x in {1..7}
do
	#calcluate sum of singleton
	> topmed_freeze2_unreklated_${x}000_singleton.txt
	for chr in {1..22}
	do
	cat temp/topmed_freeze2_${x}000.chr${chr}.overlap_removed.svm_pass.stats | awk 'BEGIN{n_sin=0} {if($3==1) n_sin=n_sin+1} END{print "chr'${chr}'\t"n_sin}' >> topmed_freeze2_unreklated_${x}000_singleton.txt
	done
	awk 'BEGIN {SUM=0} {SUM=SUM+$2} END {print "Total\t"SUM}' topmed_freeze2_unreklated_${x}000_singleton.txt >> topmed_freeze2_unreklated_${x}000_singleton.txt
done

#make sample list from bridges
# echo "sample size: 3765" > bridges_singleton.txt
# for chr in {1..22}
# do
# 	cat /net/bipolar/lockeae/final_freeze/snps/vcfs/chr${chr}/chr${chr}.filtered.sites.modified.vcf.v2.summary | grep 'PASS;AC=1-1' | awk '{print "chr'${chr}'""\t"$2}' >> bridges_singleton.txt
# done
# awk 'BEGIN {SUM=0} {SUM=SUM+$2} END {print "Total\t"SUM}' bridges_singleton.txt >> bridges_singleton.txt
