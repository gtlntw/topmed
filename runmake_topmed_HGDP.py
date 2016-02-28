#!/usr/bin/env python
#!/usr/bin/python
#!python
import sys
import os
import re
import getopt
import numpy
import subprocess

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
opts["jobNo"] = "500" #number of jobs in parallel
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
			idList = []
			infh = open(v)
			for line in infh.xreadlines():
				item = line.strip()
				if not re.search("^#+", item):
					idList.append(item)
			infh.close()
		else:
			opts["id"] = v

#clean start and end log
os.system("rm -rf log/start*.OK log/end*.OK slurm_scripts/*.*")

#if call by id list then submit up to 10 samples at a time
os.system("echo running {jobNo} parallel jobs".format(**opts))
os.system("echo submitting up to {sampleNo} individuals at a time\n".format(**opts))
try:
	idList  #check if the id list exists
	for x in range(0, len(idList), int(opts["sampleNo"])):
		list_temp = idList[x:(x+int(opts["sampleNo"]))]
		f = open("id.txt", 'w')
		f.write("# running " + str(x+1) + "th to "+ str(x+int(opts["sampleNo"])) +"th sample\n")
		f.writelines("\n".join(list_temp))
		f.close()
		os.system("echo running " + str(x+1) + "th to "+ str(x+int(opts["sampleNo"])) +"th sample")
		# print "python makefile_topmed_HGDP.py -l {launchMethod} -j {job_name} -i {id}".format(**opts)
		# print "make -f makefile_{job_name} -j {jobNo} --keep-going ".format(**opts)
		os.system("python makefile_topmed_HGDP.py -l {launchMethod} -j {job_name} -i {id}".format(**opts))
		os.system("make -f makefile_{job_name} -j {jobNo} --keep-going ".format(**opts))
		os.system("echo \n")
except NameError: ##when no id list is given
	os.system("python makefile_topmed_HGDP.py -l {launchMethod} -j {job_name} -i {id}".format(**opts))
	os.system("make -f makefile_{job_name} -j {jobNo} --keep-going ".format(**opts))


