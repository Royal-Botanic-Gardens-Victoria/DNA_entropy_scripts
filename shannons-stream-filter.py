#!/usr/bin/env python3

#theo allnutt 2021
#r3.7 calculates Shannon's information index for sample of reads, only outputs reads over threshold.
#usage
#shannons-folder.py "reads/*" filtered/ 0.8 12
#filtered/=output folder, 0.8= threshold, 12=threads
#accepts fasta, fastq and .gz
#n.b. for word size see:
#Akhter, S., Bailey, B., Salamon, P. et al. Applying Shannon's information theory to bacterial and phage genomes and metagenomes. Sci Rep 3, 1033 (2013). https://doi.org/10.1038/srep01033



import sys
import os
import re
import glob
import subprocess as sp
#import multiprocessing
import concurrent.futures
from Bio import SeqIO
import gzip
import math
import collections
from multiprocessing import Manager

digits = re.compile(r'(\d+)')
def tokenize(filename):
	return tuple(int(token) if match else token for token, match in ((fragment, digits.search(fragment)) for fragment in digits.split(filename)))

def estimate_shannon_entropy(dna_sequence):
	
	seq=dna_sequence.upper()
	if len(seq)-seq.count("N") <=0:
		a=0
		t=0
		c=0
		g=0
	else:
		a=float(seq.count("A")/(len(seq))) #-seq.count("N")
		t=float(seq.count("T")/(len(seq)))
		c=float(seq.count("C")/(len(seq)))
		g=float(seq.count("G")/(len(seq)))
	
	if a!=0:
		a1=a * math.log(a,4)
	else: 
		a1=0
	if t!=0:
		t1=t * math.log(t,4)
	else: 
		t1=0
	if c!=0:
		c1=c * math.log(c,4)
	else: 
		c1=0
	if g!=0:
		g1=g * math.log(g,4)
	else: 
		g1=0		
	
	
	entropy_value = -1 * ( a1 + t1 + c1 + g1 )
	
	#print(a,t,c,g, entropy_value)
	
	return entropy_value

def shan(i,chunk):
	
	name1=i.split("/")[-1]
	print(name1)

	f=open(i,'r')
	outfile=open(outfolder+"/"+name1,'w')
	
	n=0
	c=0
	for x in SeqIO.parse(f,'fasta'):
		
		slen=len(x.seq)
		
		#remainder=1+int(chunk*float("."+str(float(slen/chunk)).split(".")[-1]))
		goodseq=""
		
		for y in range(0,slen,chunk):
			n=n+1
			
			if y < slen-chunk:
				checkseq=str(x.seq)[y:y+256]
				shannon = estimate_shannon_entropy(checkseq)			
			
			else:
				checkseq=str(x.seq)[y:]
				shannon = estimate_shannon_entropy(checkseq)
				
			
			if shannon >=shan_limit:
				c=c+1
				goodseq=goodseq+checkseq
			
			#print(shannon,checkseq)
			
		#print(x.id,goodseq)		

		
		outfile.write(">"+str(x.description)+"\n"+str(goodseq)+"\n")
		
	print(i,"chunks=",n,"passed=",c,"%passed=",float((c/n)*100))
	
	outfile.close()

	
folder=sys.argv[1] 
filelist=glob.glob(folder)
filelist.sort(key=tokenize)
print(filelist)
outfolder = sys.argv[2] #outfolder
shan_limit=float(sys.argv[3])

chunk=int(sys.argv[4])
threads=int(sys.argv[5])

if __name__ == '__main__':


	executor = concurrent.futures.ProcessPoolExecutor(threads)
	futures = [executor.submit(shan, i,chunk) for i in filelist]
	concurrent.futures.wait(futures)

	#for i in filelist:
		#shan(i,chunk)


		