#!/bin/bash
for chr in {22..1}
do
   echo "running chromosome: ${chr}"
   nohup /net/snowwhite/home/khlin/bin/vcftools --bcf /net/topmed4/working/hmkang/freeze3a/v4/release/passgt.minDP10/ALL.chr${chr}.freeze3a.pass.gtonly.genotypes.bcf --keep /net/topmed2/working/khlin/analysis/2016-12-08_singleton_analysis_only_topmed/temp/unrelated_1.txt --remove-indels --singletons --out african.chr${chr} &
done