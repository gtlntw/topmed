#!/bin/bash

# chr=22
# id=NWD104907

for chr in `seq 1 22`;
do
	#ld-pruning 
	/net/snowwhite/home/khlin/tools/plink_1.9/plink --vcf /net/topmed2/working/khlin/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz --indep-pairwise 50 5 0.1 -out /net/topmed2/working/khlin/global_ancestry/HGDP/HGDP_chr${chr}
	#generate bim file for allele order and position
	/net/snowwhite/home/khlin/tools/plink_1.9/plink --vcf /net/topmed2/working/khlin/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz --extract /net/topmed2/working/khlin/global_ancestry/HGDP/HGDP_chr${chr}.prune.in -make-bed -out /net/topmed2/working/khlin/global_ancestry/HGDP/HGDP_chr${chr}
	#change rsid to chr:position
	cat HGDP_chr${chr}.bim | awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > HGDP_chr${chr}.bim.tmp && mv HGDP_chr${chr}.bim.tmp HGDP_chr${chr}.bim
	#generate list of snps to keep for bcftools
	cat HGDP_chr${chr}.bim | awk '{print $1"\t"$4}' > HGDP_chr${chr}.pos
	#calculate the MAF for each sub population
	for pop in `seq 1 7`;
	do
		cat HGDP_pop_id.txt | awk -F'\t' '{if($4=='${pop}') print $2" "$2}' > HGDP_pop_${pop}_id.txt
		/net/snowwhite/home/khlin/tools/plink_1.9/plink --bfile HGDP_chr${chr} --keep HGDP_pop_${pop}_id.txt --freq --out HGDP_chr${chr}_pop_${pop} --keep-allele-order
	done
done

# for chr in `seq 1 22`;
# do
# 	#extract sample vcf and convert to bfile
# 	bcftools view /net/topmed2/working/gt-release/sftp-barnes/freeze.1a/topmed.freeze1.nhlbi.791.sftp-barnes.keep.chr${chr}.gtonly.vcf.gz --types snps -M2 --exclude-uncalled -f PASS -s NWD104907 -R HGDP_chr${chr}.pos --output-type z --output-file NWD104907_filtered_chr${chr}.vcf.gz
# 	/net/snowwhite/home/khlin/tools/plink_1.9/plink --vcf NWD104907_filtered_chr${chr}.vcf.gz -make-bed --out NWD104907_filtered_chr${chr}
# 	#create variant id
# 	cat NWD104907_filtered_chr${chr}.bim | awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > NWD104907_filtered_chr${chr}.bim.tmp && mv NWD104907_filtered_chr${chr}.bim.tmp NWD104907_filtered_chr${chr}.bim
# 	#set ref allele
# 	/net/snowwhite/home/khlin/tools/plink_1.9/plink --bfile NWD104907_filtered_chr${chr} --a2-allele HGDP_chr${chr}.bim 6 2 --make-bed --out NWD104907_filtered_setref_chr${chr}
# 	#calculate the allele count for sample
# 	/net/snowwhite/home/khlin/tools/plink_1.9/plink --bfile NWD104907_filtered_setref_chr${chr} --freq counts --out NWD104907_filtered_setref_chr${chr} --keep-allele-order
# 	# /net/snowwhite/home/khlin/tools/plink_1.9/plink --bfile NWD104907_filtered_setref_chr${chr} --freq --out NWD104907_filtered_setref_chr${chr} --keep-allele-order
# done
