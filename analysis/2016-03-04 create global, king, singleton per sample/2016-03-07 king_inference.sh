#!/bin/bash

#create unrelated list based on Hyun's 4th degree inference
awk '{if($8 > 0.0884) print}' /net/topmed3/working/hmkang/freeze2/10692/aux/hgdp/topmed_freeze2_hgdp_merged.king_1_4.degree_4.kin0 > topmed_freeze2_hgdp_merged.king_1_4.degree_2.kin0
cut -f 1 topmed_freeze2_hgdp_merged.king_1_4.degree_2.kin0 |sort |uniq > related_1.tmp
cut -f 3 topmed_freeze2_hgdp_merged.king_1_4.degree_2.kin0 |sort |uniq > related_3.tmp
comm -13 related_1.tmp related_3.tmp > related_degree_2.tmp
comm -12 related_1.tmp related_3.tmp >> related_degree_2.tmp
cat related_degree_2.tmp | sort > related_degree_2.txt
comm -13 related_degree_2.txt ../../id_10597.txt > unrelated_degree_2.txt

#choose number of samples with > 0.5, done this in R

#extract sample with european ancestry > 0.5
PATH=$PATH\:/net/snowwhite/home/khlin/bin
cmd='srun --partition=topmed  --job-name=topmedsingleton' 

for chr in {1..22}
do
bcftools view -f PASS -S unrelated_degree_2_gr_0.5_eur.txt /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.genotypes.vcf.gz | bcftools query -i 'INFO/AC==1' -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' > topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.stats &
done

#calcluate sum of singleton
> topmed_freeze2_unrelated_3553_singleton.txt
for chr in {1..22}
do
cat topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.stats | awk 'BEGIN{n_sin=0} {if($3==1) n_sin=n_sin+1} END{print "chr'${chr}'\t"n_sin}' >> topmed_freeze2_unrelated_3553_singleton.txt
done
awk 'BEGIN {SUM=0} {SUM=SUM+$2} END {print "Total\t"SUM}' topmed_freeze2_unrelated_3553_singleton.txt >> topmed_freeze2_unrelated_3553_singleton.txt



# #combine into wholegenome
# file_list=''
# for chr in {1..22}
# do
# 	file_list="${file_list} topmed.europe.vcf/topmed.europe_chr${chr}.vcf.gz"
# done
# /net/snowwhite/home/khlin/tools/bcftools-1.3/bcftools concat ${file_list} -o topmed.europe.vcf/topmed.europe_genome.vcf.gz -O z

# #convert to binary bed format
# /net/snowwhite/home/khlin/tools/plink_1.9/plink --vcf topmed.europe.vcf/topmed.europe_genome.vcf.gz --thin 0.1 --make-bed --out topmed.europe.vcf/topmed.europe_genome

# #infer unrelated list up to 2 degree
# /net/snowwhite/home/khlin/tools/king/king -b topmed.europe.vcf/topmed.europe_genome.bed --unrelated --degree 2

# #extract unrelated samples and drop genotypes
# cat kingunrelated.txt | awk '{print $1}' > unrelated.txt 
# /net/snowwhite/home/khlin/tools/bcftools-1.3/bcftools view --types snps -f PASS --min-ac 1 -S unrelated.txt -G topmed.europe.vcf/topmed.europe_genome.vcf.gz -o topmed.europe.vcf/topmed.unrelated.europe_genome.sites.vcf.gz -O z

# #count singleton per sample
# /net/snowwhite/home/khlin/tools/bcftools-1.3/bcftools stats topmed.europe.vcf/topmed.unrelated.europe_genome.sites.vcf.gz > topmed.europe.vcf/topmed.unrelated.europe_genome.sites.summary


#make sample list from bridges
echo "sample size: 3765" > bridges_singleton.txt
for chr in {1..22}
do
	cat /net/bipolar/lockeae/final_freeze/snps/vcfs/chr${chr}/chr${chr}.filtered.sites.modified.vcf.v2.summary | grep 'PASS;AC=1-1' | awk '{print "chr'${chr}'""\t"$2}' >> bridges_singleton.txt
done
awk 'BEGIN {SUM=0} {SUM=SUM+$2} END {print "Total\t"SUM}' bridges_singleton.txt >> bridges_singleton.txt
