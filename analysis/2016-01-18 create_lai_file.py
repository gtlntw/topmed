#!/usr/bin/env python
#!/usr/bin/python
#!python
import sys
import os
import re
import getopt
import numpy
import subprocess
import gzip
from collections import Counter

def usage():
	print """
	submit the pipeline by batch of samples

	 makefile [-h] [-i <input filename>] [-j <job name>] [-l <launch method>]

	 -h	  ask help

	 -i   input filename

	 -l	  launch by local or slurm

	 -j	  job name

	 """

#option variables
opts = {}
opts["id"] = "." #default is to read from id.txt
opts["launchMethod"] = "slurm"
opts["job_name"] = "topmedlai"
opts["jobNo"] = "200" #number of jobs in parallel
opts["sampleNo"] = "100" #up to number of samples to submit at a time

o, a = getopt.getopt(sys.argv[1:], "i:l:j:h")
for k,v in o:
	opts[k] = v
	if opts.has_key("-h"):
		usage(); sys.exit(0)
	elif k == "-l":
		opts["launchMethod"] = v
		if opts["launchMethod"] != "local" and opts["launchMethod"] != "slurm":
			print "Launch method has to be local or slurm"; exit(1)
	elif k == "-j":
		opts["jobName"] = v
		opts["makeFile"] = "makefile_{jobName}".format(**opts)
	elif k == "-i":
		if re.search("txt", v): #if it is a file, read the id list
			##readin id
			idList = []
			infh = open(v)
			for line in infh.xreadlines():
				item = line.strip()
				if not re.search("^#+", item):
					idList.append(item)
			infh.close()
		else:
			opts["id"] = v

#function to calculate the proportion
def laiProp(laiListIdv):
	nMarker = float(len(laiListIdv))
	ancestry = [1,2,3,4,5,6,7]
	ancestryCount = [0,0,0,0,0,0,0]
	for i in ancestry:
		ancestryCount[i-1] = laiListIdv.count(str(i))
	return [str(x/nMarker) for x in ancestryCount]

##create a list to store result
propList = []

##readin lai result and loop over the individual's chromosomes to calculate proportion
laiList = []
#idList = ['NWD100295', 'NWD100505']
for a in idList:
	opts['id'] = a
	print "Processing {id}".format(**opts)
	laiListIdv = []
	for i in xrange(1,23):
		opts['chr'] = i
		#print "Chromosome {chr}".format(**opts)
		infh = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt".format(**opts))
		for line in infh.xreadlines():
			laiPos = line.strip().split()
			[laiListIdv.append(x) for x in laiPos] #append the lai at each position
	laiList.append([a] + laiProp(laiListIdv)) #save id and proportion

##write the result into a file
f = open("/net/topmed2/working/khlin/analysis/laiProp.txt", 'w')
#print ["\t".join(x) for x in laiList]
f.writelines("id\tpop1\tpop2\tpop3\tpop4\tpop5\tpop6\tpop7\n")
f.writelines("\n".join(["\t".join(x) for x in laiList]))

