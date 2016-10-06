#!/usr/bin/env python
#!/usr/bin/python
#!python
import sys
import os
import re
import getopt
import numpy
from scipy.stats import expon
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
opts["gzip"] = "FALSE"

o, a = getopt.getopt(sys.argv[1:], "i:l:j:h:c:z:")
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
	elif k == "-z":
		opts["gzip"] = "TRUE"
	elif k == "-i":
		if re.search("txt|ids", v): #if it is a file, read the id list
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

#clean start and end log

#if call by id list then submit up to 10 samples at a time
# os.system("echo running {jobNo} parallel jobs".format(**opts))
# os.system("echo submitting up to {sampleNo} individuals at a time\n".format(**opts))
# try:
# 	idList  #check if the id list exists
# 	for x in range(0, len(idList), int(opts["sampleNo"])):
# 		list_temp = idList[x:(x+int(opts["sampleNo"]))]
# 		f = open("id.txt", 'w')
# 		f.writelines("\n".join(list_temp))
# 		f.close()
# 		os.system("echo running " + str(x+1) + "th to "+ str(x+int(opts["sampleNo"])) +"th sample")
# 		# print "python makefile_topmed_HGDP.py -l {launchMethod} -j {job_name} -i {id}".format(**opts)
# 		# print "make -f makefile_{job_name} -j {jobNo} --keep-going ".format(**opts)
# 		os.system("python makefile_topmed_HGDP.py -l {launchMethod} -j {job_name} -i {id}".format(**opts))
# 		os.system("make -f makefile_{job_name} -j {jobNo} --keep-going ".format(**opts))
# 		os.system("echo \n")
# except NameError: ##when no id list is given
# 	os.system("python makefile_topmed_HGDP.py -l {launchMethod} -j {job_name} -i {id}".format(**opts))
# 	os.system("make -f makefile_{job_name} -j {jobNo} --keep-going ".format(**opts))

def convertLAI_2_GT(x):
	lookupDict = {
		'1/1': 1,
		'1/2': 2,
		'1/3': 3,
		'1/4': 4,
		'1/5': 5,
		'1/6': 6,
		'1/7': 7,
		'2/2': 8,
		'2/3': 9,
		'2/4': 10,
		'2/5': 11,
		'2/6': 12,
		'2/7': 13,
		'3/3': 14,
		'3/4': 15,
		'3/5': 16,
		'3/6': 17,
		'3/7': 18,
		'4/4': 19,
		'4/5': 20,
		'4/6': 21,
		'4/7': 22,
		'5/5': 23,
		'5/6': 24,
		'5/7': 25,
		'6/6': 26,
		'6/7': 27,
		'7/7': 28}
	match_gt = lookupDict.get(x, -9)
	if match_gt == -9:
		x_swap = x[2] + '/' + x[0]
		match_gt = lookupDict.get(x_swap, 0)
	result = ['0' for x in range(28)]
	result[match_gt-1] = '1' #gt matches which combination
	return result

def convertGT_2_LAI(x):
	x = [float(i) for i in x]
	lookupDict = {
		0: '1/1',
		1: '1/2',
		2:'1/3',
		3:'1/4',
		4:'1/5',
		5:'1/6',
		6:'1/7',
		7:'2/2',
		8:'2/3',
		9:'2/4',
		10:'2/5',
		11:'2/6',
		12:'2/7',
		13:'3/3',
		14:'3/4',
		15:'3/5',
		16:'3/6',
		17:'3/7',
		18:'4/4',
		19:'4/5',
		20:'4/6',
		21:'4/7',
		22:'5/5',
		23:'5/6',
		24:'5/7',
		25:'6/6',
		26:'6/7',
		27:'7/7'}
	return lookupDict.get(x.index(max(x)))
	
# def convertLAI_2_GT_LI(lai_left, lai_right, w):
# 	if(lai_left == lai_right):
# 		return convertLAI_2_GT(lai_left)
# 	left_result = convertLAI_2_GT(lai_left)
# 	right_result = convertLAI_2_GT(lai_right)
# 	result = [str(w*float(right_result[x]) + (1-w)*float(left_result[x])) for x in xrange(28)]
# 	return result

def convertGT_LI(lai_left, lai_right, w):
	if(lai_left == lai_right):
		return lai_left
	result = [str(w*float(lai_right[x]) + (1-w)*float(lai_left[x])) for x in xrange(28)]
	return result

#readin global ancestry
print "Reading in Global ancestry\n".format(**opts)
gaiList = []
infh = open("/net/topmed2/working/khlin/global.txt".format(**opts))
for line in infh.xreadlines():
	item = line.strip().split(' ')
	if item[0] in idList: #assume id.txt and global.txt are both sorted
		gaiList.append(item[1:8])
infh.close()

#convert global ancestry to dosage
gaiListGT = [[str(float(gaiList[z][x])*float(gaiList[z][y])) for x in xrange(0,7) for y in xrange(x,7)] for z in xrange(len(gaiList))]

##readin snp info
print "processing chromosome {chr}".format(**opts)

##readin lai result
laiList = []
print "Reading in LAI results".format(**opts)
for i in idList:
	opts['id'] = i
	laiListIdv = []
	posListHgdp = []
	posCommon = []
	infh = open("/net/topmed2/working/khlin/common_site/topmed_freeze12_chr{chr}_HGDP_common.txt".format(**opts))
	for line in infh.xreadlines():
		item = line.strip().split('\t')
		posCommon.append(item[1])
	infh.close()
	infh = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.pos".format(**opts))
	for line in infh.xreadlines():
		item = line.strip().split('\t')
		posListHgdp.append(item[1])
	infh.close()
	infh = open("/net/topmed2/working/khlin/output/LAI/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt".format(**opts))
	for a, line in enumerate(infh.xreadlines()):
		if posListHgdp[a] in posCommon:
			temp = line.strip()
			item = '/'.join(temp.split())
			if int(item[0]) > int(item[2]):
				item = item[2]+'/'+item[0] #e.g. change 4/1 to 1/4 
			laiListIdv.append(item)
	laiList.append(laiListIdv)
	infh.close()

##convert local ancestry result to dosage
laiListGT = [[convertLAI_2_GT(laiList[i][j]) for j in xrange(len(laiList[0]))] for i in xrange(len(laiList))]

##output vcf file
if opts['gzip'] == "TRUE": f = gzip.open("/net/topmed2/working/khlin/vcf/chr{chr}.lai.vcf.gz".format(**opts), 'w', compresslevel=1)
else: f = open("/net/topmed2/working/khlin/vcf/chr{chr}.lai.vcf".format(**opts), 'w')
f.writelines("##fileformat=VCFv4.1\n")
f.writelines("##contig=<ID={chr}>\n".format(**opts))
f.writelines("##FILTER=<ID=LI,Description=\"Linear interpolation for SNPs lying between two HGDP markers\"\n")
f.writelines("##FILTER=<ID=EDGE,Description=\"SNPs before/after the first/last HGDP marker. Weighted mean of global ancesrtry and the local ancestry of the closest HGDP marker.\"\n")
f.writelines("##FILTER=<ID=HGDP,Description=\"HGDP markers\"\n")
f.writelines("##FORMAT=<ID=COMB,Number=1,Type=String,Description=\"Local Ancestry combination. 1:Sub-Saharan Africa, 2:Central and South Asia, 3:East Asia, 4:Europe, 5:Native America, 6:Oceania, 7:West Asia and North Africa.\">\n".format(**opts))
f.writelines("##FORMAT=<ID=COMB_DOSAGE,Number=28,Type=Float,Description=\"The likelihood of Local Ancestry combinations appearing in the order of 1/1,1/2,1/3,1/4,1/5,1/6,1/7,2/2,2/3,2/4,2/5,2/6,2/7,3/3,3/4,3/5,3/6,3/7,4/4,4/5,4/6,4/7,5/5,5/6,5/7,6/6,6/7,7/7\">\n".format(**opts))
# f.writelines("##FORMAT=<ID=1/1,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=1/2,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=1/3,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=1/4,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=1/5,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=1/6,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=1/7,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=2/2,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=2/3,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=2/4,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=2/5,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=2/6,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=2/7,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=3/3,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=3/4,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=3/5,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=3/6,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=3/7,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=4/4,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=4/5,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=4/6,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=4/7,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=5/5,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=5/6,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=5/7,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=6/6,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=6/7,Number=1,Type=Float>,Description=\"Dosage\"\n")
# f.writelines("##FORMAT=<ID=7/7,Number=1,Type=Float>,Description=\"Dosage\"\n")

#input the position of topmed snps
print "Reading in Topmed SNP position".format(**opts)
posListTopmed = []
infh = open("/net/topmed2/working/khlin/vcf/chr{chr}.freeze2.pos".format(**opts))
for line in infh.xreadlines():
	item = line.strip().split('\t')
	posListTopmed.append(item)
infh.close()

#add a decsription of genotype allele
hgdp_snp_idx = 0
hgdp_snp_pos = int(posCommon[hgdp_snp_idx])
n_sample = len(laiList)
n_topmed_snp = len(posListTopmed)
first_hgdp_snp_pos = int(posCommon[0])
last_hgdp_snp_pos = int(posCommon[-1])
f.writelines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(idList) + "\n")
for h in xrange(n_topmed_snp):
	if h % 100000 == 0: print str(h)+"/"+str(n_topmed_snp)+" SNPs done"
	topmed_snp_pos = int(posListTopmed[h][1])
	if(topmed_snp_pos < first_hgdp_snp_pos):
		gtList = []
		#need to employ weight to the closest HGDP marker
		w = 1 - expon.cdf(first_hgdp_snp_pos - topmed_snp_pos, scale=5000000)  # p(recom=FALSE) mean length is 5MB
		for i in xrange(n_sample):
			# print str(i) + " " + str(h)
			# if (first_hgdp_snp_pos - topmed_snp_pos) > 5000000: #use 5MB as cut off
			# 	result = gaiListGT[i]
			# else:
			# 	result = laiListGT[i][0]
			result = convertGT_LI(gaiListGT[i], laiListGT[i][0], w)
			# result = [gaiListGT[i][j]*w + (1-w)*laiListGT[i][0] for j in xrange(7)] 
			gtList.append(convertGT_2_LAI(result)+":"+','.join(result))
		f.writelines("\t".join(posListTopmed[h]) + "\t.\t.\t.\t.\tEDGE\t.\tCOMB:COMB_DOSAGE\t" + "\t".join(gtList) + "\n") #
		continue
	if(topmed_snp_pos > last_hgdp_snp_pos):
		gtList = []
		#need to employ weight to the closest HGDP marker
		w = 1 - expon.cdf(topmed_snp_pos - last_hgdp_snp_pos, scale=5000000)  # p(recom=FALSE) mean length is 5MB
		for i in xrange(n_sample):
			# print str(i) + " " + str(h)
			# if (topmed_snp_pos - last_hgdp_snp_pos) > 5000000: #use 5MB as cut off
			# 	result = gaiListGT[i]
			# else:
			# 	result = laiListGT[i][-1]
			result = convertGT_LI(gaiListGT[i], laiListGT[i][-1], w)
			gtList.append(convertGT_2_LAI(result)+":"+','.join(result))
		f.writelines("\t".join(posListTopmed[h]) + "\t.\t.\t.\t.\tEDGE\t.\tCOMB:COMB_DOSAGE\t" + "\t".join(gtList) + "\n") #
		continue
	if(topmed_snp_pos == hgdp_snp_pos):
		#print "matched HGDP position"
		gtList = []
		for i in xrange(n_sample):
			# print str(i) + " " + str(h)
			gtList.append(laiList[i][hgdp_snp_idx]+":"+','.join(laiListGT[i][hgdp_snp_idx]))
		f.writelines("\t".join(posListTopmed[h]) + "\t.\t.\t.\t.\tHGDP\t.\tCOMB:COMB_DOSAGE\t" + "\t".join(gtList) + "\n") #
		hgdp_snp_idx += 1 #move to the next hgdP snp
		try: #prevent out of range for the last snp
			hgdp_snp_pos = int(posCommon[hgdp_snp_idx])
		except IndexError:
			print "reach the end"
			continue
	else:
		#print "matched HGDP position"
		#!check the position residing between two hgdp snps
		gtList = []
		w = (float(topmed_snp_pos)-float(posCommon[hgdp_snp_idx-1]))/(float(posCommon[hgdp_snp_idx]) - float(posCommon[hgdp_snp_idx-1])) #weight by the proximity to the right flanking snp
		for i in xrange(n_sample):
			# print str(i) + " " + str(h)
			result = convertGT_LI(laiListGT[i][hgdp_snp_idx-1], laiListGT[i][hgdp_snp_idx], w)
			gtList.append(convertGT_2_LAI(result)+":"+','.join(result))
		f.writelines("\t".join(posListTopmed[h]) + "\t.\t.\t.\t.\tLI\t.\tCOMB:COMB_DOSAGE\t" + "\t".join(gtList) + "\n") #
f.close()
print "\n"
