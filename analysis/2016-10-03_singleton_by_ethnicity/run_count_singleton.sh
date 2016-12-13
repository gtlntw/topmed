#!/bin/bash
#create 2nd unrelated list
# king -b /net/topmed4/working/hmkang/freeze3a/v4/hgdp/topmed_freeze3a_hgdp_autosomes_20160816.plink.bed --unrelated --degree 2

#create the global ancestry grouping based on >0.5 and max rules
python calculate_ethnicity.py

#create makefile to count singleton by ethnicity, increasing size, chr
python calculate_singleton_by_ethnicity.py

#run the pipeline
nohup make -j 30 > log.txt 2>&1 &
