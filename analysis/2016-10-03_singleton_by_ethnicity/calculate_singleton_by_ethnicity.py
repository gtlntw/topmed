#!/usr/bin/env python
#!/usr/bin/python
#!python

import sys
import os
import re
import getopt
import numpy
from string import Template

######################
#Initialization
######################
idList=[]
unrelatedList=[]
ethnicityList=[] #unrelated samples only

def find_gt_half(list):
	idx = 0
	for i in xrange(len(list)):
		if float(list[i]) > 0.5:
			idx = i + 1 #save as 1,2,3,4,5,6,7
	return idx

######################
#read in ethnicity ancestry and unrelated list
######################
infh = open("unrelated_degree_2.txt")
for line in infh.xreadlines():
	item = line.strip()
	unrelatedList.append(item)
infh.close()

infh = open("/net/topmed2/working/khlin/ethnicity.txt")
for line in infh.xreadlines():
	sample_id, ethnicity = line.strip().split("\t")
	if sample_id in unrelatedList: #only unrelated
		ethnicityList.append([sample_id, ethnicity])
infh.close()
print ethnicityList
print len(ethnicityList)

######################
#output unrelated sample list by ancestry
######################
import random
#set seed, make sure results are tractable
random.seed(464351)
#shuffle the ethnicity list
random.shuffle(ethnicityList)
#output each ethnicity
for i in xrange(8):
	outfh = open("temp/unrelated_{ethnicity}.txt".format(ethnicity=i), 'w')
	for j in ethnicityList:
		if int(j[1]) == i:
			outfh.writelines("\t".join(j)+"\n")
	outfh.close()

######################
#count singleton
######################
import subprocess
resultList=[]
idList=[]
size=0
chr=22
ethnicityIdx=1

infh = open("temp/unrelated_{ethnicityIdx}.txt".format(ethnicityIdx=ethnicityIdx))
for line in infh:
	items = line.strip().split()
	idList.append(items[0])
infh.close()
outfh = open("temp/id_temp.txt", 'w')
outfh.writelines('\n'.join(idList[0:1000]))
outfh.close()

cmd = "bcftools view --type snps -G -S temp/id_temp.txt /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr{chr}.overlap_removed.svm_pass.genotypes.vcf.gz | bcftools query -i 'INFO/AC==1' -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' | wc -l".format(chr=chr)
# os.system(cmd)
result = subprocess.check_output(cmd, shell=True)
print result