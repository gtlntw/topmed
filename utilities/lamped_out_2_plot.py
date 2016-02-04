#!/usr/bin/env python
#!/usr/bin/python
#!python

import sys
import getopt

def usage():
	print """
	convert lampped out format to one haplotype one column per position

	 lamped_out_2_plot [-h] [-i <input filename>] [-o <output filename>]

	 -h                  print this message
	 
	 -i                  input lamped filename

	 -p                  position filename
	        
	 -o                  output filenam

	 """

o, a = getopt.getopt(sys.argv[1:], 'i:o:p:h')
opts = {}
for k,v in o:
	opts[k] = v
	if opts.has_key('-h'):
	    usage(); sys.exit(0)
	if k == "-i":
		inputName = v
	if k == "-o":
		outputName = v	
	if k == "-p":
		posName = v	

##Read in lamped out file
infh= open(inputName)
for line in infh.xreadlines():
    fields = line.split()
dosage=[]
switch_pos=[]
for item in fields:
	a, b=((item.split(':')))
	dosage.append(a)
	switch_pos.append(b)
infh.close()
# print(switch_pos)
##Read in position file
infh= open(posName)
pos=[]
for line in infh.xreadlines():
    pos.append(line.split()[1])
infh.close()
#dictionary convert output 1,2,3,4,5,6,7 to diploid format 1 1, 1 2 ...... 1: african  4: europe  5: native american
mapping = {'00': '4 4', '01': '4 5', '02': '4 1', '11': '5 5', '12': '5 1', '22': '1 1'}

##for each position write out the local ancestry
ofh = open(outputName, "w")
start_pos_index = 0
for idx, value in enumerate(switch_pos):
	switch_pos_index =int(value)-1
	ancestry = mapping[str(dosage[idx])]
	for i in pos[start_pos_index:switch_pos_index]:
		ofh.write(ancestry+"\n")
	start_pos_index = switch_pos_index;
ofh.write(mapping[str(dosage[-1])]+"\n") #print the last snp
ofh.close()
