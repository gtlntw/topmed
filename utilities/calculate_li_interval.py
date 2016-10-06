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


#calculate the number of LI interval
n = len(idList)
meanSwitchIntervalList = []
snpLiIntervalList = []
meanSwitchSnpList = []
chrList = []
for chr in xrange(1,23):
	print "processing chromosome: " + str(chr)
	opts['chr'] = chr
	posCommon = []
	infh = open("/net/topmed2/working/khlin/common_site/topmed_freeze12_chr{chr}_HGDP_common.txt".format(**opts))
	for line in infh.xreadlines():
		item = line.strip().split('\t')
		posCommon.append(float(item[1]))
	infh.close()
	print len(posCommon)

	#calculate the number of snp in each LI interval
	countLiInterval = []
	hgdpIdx = 1
	countTopmedSnp = 0
	infh = open("/net/topmed2/working/khlin/vcf/chr{chr}.freeze2.pos".format(**opts))
	for line in infh.xreadlines():
		item = float(line.strip().split('\t')[1])
		if posCommon[hgdpIdx] < item: #move to the next interval
			countLiInterval.append(countTopmedSnp)
			countTopmedSnp = 0
			hgdpIdx +=1
		try:
			if posCommon[hgdpIdx-1] < item and item < posCommon[hgdpIdx]:
				countTopmedSnp +=1
		except IndexError:
			print "reach the end"
			break
	infh.close()
	print len(countLiInterval)

	#run by individual, read in LAI and determine whether there is a switch, save
	posListHgdp = []
	countSwitchInterval = [0.] * (len(posCommon)-1)
	for i in idList:
		print "processing " + i
		opts['id'] = i
		infh = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.pos".format(**opts))
		for line in infh.xreadlines():
			item = line.strip().split('\t')
			posListHgdp.append(float(item[1]))
		infh.close()
		infh = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt".format(**opts))
		commonIdx = 0 # start from the first common snp
		for idx, line in enumerate(infh.xreadlines()):
			try:
				posCommon[commonIdx]
			except IndexError:
				print "reach the end"
				commonIdx -= 1	
			if posListHgdp[idx] == posCommon[commonIdx]:
				commonIdx +=1
				temp = line.strip()
				item = '/'.join(temp.split())
				if float(item[0]) > float(item[2]):
					item = item[2]+'/'+item[0] #e.g. change 4/1 to 1/4
				if idx == 0: #initialize the first
					preLai = item
				if item != preLai:
					print commonIdx-1
					countSwitchInterval[commonIdx-2] += 1
					preLai = item
		infh.close()

	# print len(idList) #number of individuals\
	chrList.append(chr)
	snpLiIntervalList.append(sum(countLiInterval)) #number of snps in intervals
	meanSwitchIntervalList.append(sum(countSwitchInterval)/n) #mean number of switch
	meanSwitchSnpList.append(sum([a*b for a,b in zip(countLiInterval, countSwitchInterval)])/n) #mean LI snps per individuals

print snpLiIntervalList
print meanSwitchIntervalList
print meanSwitchSnpList

outf = open("calculate_li_interval.csv", 'w')
outf.writelines("chr,n,snpLiInterval,meanSwitchInterval,meanSwitchSnp\n")
for a,b,c,d in zip(chrList, snpLiIntervalList, meanSwitchIntervalList, meanSwitchSnpList):
	outf.writelines(",".join([str(a),str(n),str(b),str(c),str(d)])+"\n")
outf.close()
##readin snp info
# distList = []
# for chrom in xrange(22,23):
# 	opts['chr'] = chrom
# 	print "processing chromosome {chr}".format(**opts)

# 	##readin lai result
# 	laiList = []
# 	for i in idList:
# 		opts['id'] = i
# 		print "Processing {id}".format(**opts)
# 		# print "Reading in HGDP SNP position".format(**opts)
# 		posListHgdp = []
# 		infh = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.pos".format(**opts))
# 		for line in infh.xreadlines():
# 			item = line.strip().split('\t')
# 			posListHgdp.append(item[1])
# 		infh.close()
		
# 		# print "Reading in LAI resuls".format(**opts)
# 		laiListIdv = []
# 		infh2 = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt".format(**opts))
# 		for line in infh2.xreadlines():
# 			item = line.strip().split(' ')
# 			laiListIdv.append(item)
# 		infh2.close()	

# 		lai_current = laiListIdv[0]
# 		lai_left_position = float(posListHgdp[0])
# 		lai_right_position = float(posListHgdp[0])
# 		dist = lai_right_position - lai_left_position
# 		for i, ele in enumerate(laiListIdv):
# 			if not ((lai_current[0] == ele[0] and lai_current[1] == ele[1]) or (lai_current[0] == ele[1] and lai_current[1] == ele[0])): #if switch happen
# 				lai_current = ele
# 				lai_right_position = float(posListHgdp[i-1])
# 				dist = lai_right_position - lai_left_position
# 				# if dist != 0: distList.append(dist)
# 				distList.append(dist)
# 				lai_left_position = float(posListHgdp[i])
# 		##Handle the last segment
# 		lai_right_position = float(posListHgdp[-1])
# 		dist = lai_right_position - lai_left_position
# 		# if dist != 0: distList.append(dist)
# 		distList.append(dist)

#         # print str(i) + " " + str(len(laiListIdv))
# print "mean distance is {:,f}".format(numpy.mean(distList))
