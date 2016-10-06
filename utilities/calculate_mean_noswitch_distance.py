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

if 'idList' not in globals():
	sys.exit("need to input id list: -i xxxx.txt")

##readin snp info
distList = []
for chrom in xrange(1,23):
	opts['chr'] = chrom
	print "processing chromosome {chr}".format(**opts)

	##readin lai result
	laiList = []
	for i in idList:
		opts['id'] = i
		print "Processing {id}".format(**opts)
		# print "Reading in HGDP SNP position".format(**opts)
		posListHgdp = []
		infh = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.pos".format(**opts))
		for line in infh.xreadlines():
			item = line.strip().split('\t')
			posListHgdp.append(item[1])
		infh.close()
		
		# print "Reading in LAI resuls".format(**opts)
		laiListIdv = []
		infh2 = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt".format(**opts))
		for line in infh2.xreadlines():
			item = line.strip().split(' ')
			laiListIdv.append(item)
		infh2.close()	

		lai_current = laiListIdv[0]
		lai_left_position = int(posListHgdp[0])
		lai_right_position = int(posListHgdp[0])
		for i, ele in enumerate(laiListIdv):
			if not ((lai_current[0] == ele[0] and lai_current[1] == ele[1]) or (lai_current[0] == ele[1] and lai_current[1] == ele[0])): #if switch happen
				lai_current = ele
				lai_right_position = int(posListHgdp[i-1])
				dist = lai_right_position - lai_left_position
				# if dist != 0: distList.append(dist)
				distList.append(dist)
				lai_left_position = int(posListHgdp[i])
		##Handle the last segment
		lai_right_position = int(posListHgdp[-1])
		dist = lai_right_position - lai_left_position
		# if dist != 0: distList.append(dist)
		distList.append(dist)

        # print str(i) + " " + str(len(laiListIdv))
print "mean distance is {:,f}".format(numpy.mean(distList))
