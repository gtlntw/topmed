#combine vcf
bcftools merge /net/topmed3/working/hmkang/freeze2/10597.v2/subset/sftp-ellinor/topmed_freeze2.chr20.svm_pass_gt.sftp-ellinor.vcf.gz /net/topmed3/working/hmkang/freeze2/10597.v2/subset/sftp-ramachandran/topmed_freeze2.chr20.svm_pass_gt.sftp-ramachandran.vcf.gz -o ellinor_ramachandran_chr20.vcf.gz -O z

#convert to binary bed format
/net/snowwhite/home/khlin/tools/plink_1.9/plink --vcf ellinor_ramachandran_chr20.vcf.gz --make-bed --out ellinor_ramachandran_chr20

#infer unrelated list up to 2 degree
/net/snowwhite/home/khlin/tools/king/king -b ellinor_ramachandran_chr20.bed --unrelated --degree 2

#subset the sample
 cut -f 1 kingunrelated.txt > ellinor_ramachandran_chr20_subset_id.txt
bcftools view -S ellinor_ramachandran_chr20_subset_id.txt ellinor_ramachandran_chr20.vcf.gz -O b | bcftools stats > ellinor_ramachandran_chr20_subset_stats.txt


#singleton count


#make sample list from bridges
bcftools query -l /net/bipolar/lockeae/final_freeze/snps/vcfs/chr20/chr20.filtered.20150702.KEEP.vcf.gz > bridges_id.txt
#then I shuffle the sample list in Sublime Text
#extract subsample from bridges ps:bcftools gave me trouble!! so I used vcftools
vcftools --gzvcf /net/bipolar/lockeae/final_freeze/snps/vcfs/chr20/chr20.filtered.20150702.KEEP.vcf.gz --keep bridges_id.txt --recode --stdout | gzip -c > bridges.subset.chr20.vcf.gz
#singleton count
bcftools stats /net/bipolar/lockeae/final_freeze/snps/vcfs/chr20/chr20.filtered.20150702.KEEP.vcf.gz > bridges_stats.txt
