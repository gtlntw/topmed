#!/usr/bin/env python
#!/usr/bin/python
#!python
import sys
import getopt

def usage():
	print """
	convert lampped out format to one haplotype one column per position

	 lamped_out_2_plot [-h] [-i <input filename>] [-p <position filename>] [-o <output filename>]

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


#dictionary convert output 1: euro 2: africa 3: native, to diploid format 1 1, 1 2 ...... 1: african  4: europe  5: native american
mapping = {'11': '4 4', '12': '4 1', '13': '4 5', '22': '1 1', '23': '1 5', '33': '5 5'}

##Read in eila out file
infh= open(inputName)
ofh = open(outputName, "w")
for line in infh.xreadlines():
	if(line.strip() != 'NA'): 
		ancestry = mapping[str(line.strip())]
		ofh.write(ancestry+"\n")
	else:
		ofh.write(ancestry+"\n") #when NA just output the ancestry at the last position
ofh.close()