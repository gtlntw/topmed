#!/usr/bin/env python
#!/usr/bin/python
#!python

import sys
import os
import re
import getopt
import numpy
import subprocess
import collections

singleton_indv = []

######################
#read in singleton count
######################
##usingleton count
for chr in xrange(1,23):
	infh = open("african.chr{chr}.singletons".format(chr=chr))
	next(infh)
	for line in infh.xreadlines():
		item = line.strip().split()
		singleton_indv.append(item[4])
	infh.close()

counter=collections.Counter(singleton_indv)
print(counter)

######################
#output the count of singleton for each individual
######################
##usingleton count
outfh = open("count.txt", 'w')
for x in counter:
	output = x + "\t" + str(counter[x]) + "\n"
	outfh.writelines(output) #first i samples
outfh.close()