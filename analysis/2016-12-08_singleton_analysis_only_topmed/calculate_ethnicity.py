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

def find_gt_half(list):
	idx = 0
	for i in xrange(len(list)):
		if float(list[i]) > 0.5:
			idx = i + 1 #save as 1,2,3,4,5,6,7
	return idx

######################
#read in global ancestry
######################
infh = open("global_freeze3a_topmed_only.txt")
outfh = open("ethnicity_half_topmed_only.txt", 'w')
outfh1 = open("ethnicity_max_topmed_only.txt", 'w')
for line in infh.xreadlines():
	items = line.strip().split(" ")
	id = items[0]
	ancestry_prob = [float(x) for x in items[1:8]]
	#find >0.5 ancestry
	ethnicity_half = find_gt_half(ancestry_prob)
	# print ethnicity
	outfh.writelines([str(items[0])+"\t"+str(ethnicity_half)+"\n"])

	#find max ancestry
	print ancestry_prob
	max_value = max(ancestry_prob)
	ethnicity_max = ancestry_prob.index(max_value) + 1 #save as 1,2,3,4,5,6,7
	# print ethnicity_max
	outfh1.writelines([str(items[0])+"\t"+str(ethnicity_max)+"\n"])
