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
opts["chr"] = "22"

o, a = getopt.getopt(sys.argv[1:], "i:l:j:h:c:")
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
	elif k == "-c":
		opts["chr"] = v
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


##readin freeze 1
for chrom in xrange(1,23):
	opts['chr'] = chrom
	freeze1 = []
	infh = open("topmed_freeze1_chr{chr}_HGDP_common.txt".format(**opts))
	for line in infh.xreadlines():
		item = line.strip().split('\t')
		freeze1.append(int(item[1]))
	infh.close()

	freeze2 = []
	infh = open("topmed_freeze2_chr{chr}_HGDP_common.txt".format(**opts))
	for line in infh.xreadlines():
		item = line.strip().split('\t')
		freeze2.append(int(item[1]))
	infh.close()

	freeze12 = sorted(list(set.intersection(set(freeze1), set(freeze2))))

	f = open("topmed_freeze12_chr{chr}_HGDP_common.txt".format(**opts), 'w')
	for i in freeze12:
		f.writelines("{chr}\t".format(**opts)+str(i)+"\n")
	f.close()	