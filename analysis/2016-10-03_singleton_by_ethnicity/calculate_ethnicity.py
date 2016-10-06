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
infh = open("/net/topmed2/working/khlin/global.txt")
outfh = open("/net/topmed2/working/khlin/ethnicity.txt", 'w')
for line in infh.xreadlines():
	items = line.strip().split(" ")
	ethnicity = find_gt_half(items[1:7])
	# print ethnicity
	outfh.writelines([str(items[0])+"\t"+str(ethnicity)+"\n"])