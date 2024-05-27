#!/usr/bin/env python3

#theo allnutt 2021
#see r10.2 for notes
#outputs Shannon's and kmer count entropy measures for DNA sequence file or files
#usage:
#kmercount-shannons.py "input-folder/*" output-folder/ <kmerlen> <threads>
#kmercount-shannons.py "mega353.fasta" ./ 8 1
#kmercount-shannons.py "reads/*.fasta" kmercounts/ 8 12
#kmercount-shannons.py "reads/*.fastq.gz" kmercounts/ 4 12


import sys
import os
import re
import glob
import subprocess as sp
import concurrent.futures
from Bio import SeqIO
import gzip
import math
import collections


def tokenize(filename):
	digits = re.compile(r'(\d+)')
	return tuple(int(token) if match else token for token, match in ((fragment, digits.search(fragment)) for fragment in digits.split(filename)))

def estimate_shannon_entropy(dna_sequence):
	
	seq=dna_sequence.upper()
	if len(seq)-seq.count("N") <=0:
		a=0
		t=0
		c=0
		g=0
	else:
		a=float(seq.count("A")/(len(seq)-seq.count("N")))
		t=float(seq.count("T")/(len(seq)-seq.count("N")))
		c=float(seq.count("C")/(len(seq)-seq.count("N")))
		g=float(seq.count("G")/(len(seq)-seq.count("N")))
	
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

def kmercount(dna_sequence,wordsize):
	
	wordsizeint=int(wordsize)
	kmers=[]
	seqlen=float(len(dna_sequence))
	wordsize=float(wordsize)
	
	
	#n.b. this posscount is modelled for shorter sequences where b^len count is not true
	posscount=float((4**wordsize)*(1-math.exp(-(1/(4**wordsize))*seqlen)))
	
	#use infinate len seq posscount
	#posscount=4**wordsize
	
	#print(posscount)
	for x in range(len(dna_sequence)-wordsizeint):
		
		kbit=dna_sequence[x:x+wordsizeint]
		if "N" not in kbit:
			
		#if kbit not in kmers:
			#print(x,kbit)
			kmers.append(kbit) #also try if to test before set see which is faster: 
		
		#if len(kmers)>=posscount:
			#break
	kmers=set(kmers)
	
	kratio=float(len(kmers)/posscount)
	#kratio is ratio of kmer count to maximum possilbe count for the sequence length
	return (len(kmers),kratio)
		

def shan(i,wordn):
	
	name0=i.split("/")[-2]
	name1=i.split("/")[-1]
	name1=name0+"_"+name1
	#print(name1)
	ext1=name1.split(".")[-1]
	
	if ext1=="gz":
		#print("gz detected")
		f=gzip.open(i,'rt')
		k=name1.split(".")[-2]
		name1=name1[:-3]
	else:
		f=open(i,'r')
		k=ext1
	
	if k[-1]=="a":
		fmt="fasta"
		#print("fasta")
	if k[-1]=="q":
		fmt="fastq"
		#print("fastq")
	c=0
	
	st=0
	outfile=open(outfolder+"/"+name1+".shan",'w')
	outfile.write("name\tbp\tShannon's\tkmer_count\tkmer_ratio\n")
	
	for x in SeqIO.parse(f,fmt):
		c=c+1
		
		if c<=sample:
		
			shannon = estimate_shannon_entropy(str(x.seq))			
			kcount=kmercount(str(x.seq),wordn)
			
			outfile.write(str(x.description)+"\t"+str(len(x.seq))+"\t"+str(shannon)+"\t"+str(kcount[0])+"\t"+str(kcount[1])+"\n")
			
			#print(x.id,'len=',len(x.seq),'S=',shannon,'nK=',kcount[0],'rK=',kcount[1])
			st=st+shannon
		else:
			break
		
	shanmean=float(st/c)
	
	print(i,c,"reads",shanmean,"mean Shannons")
	
	outfile.close()

#global

folder=sys.argv[1] 
filelist=glob.glob(folder)
filelist.sort(key=tokenize)
print(filelist)
outfolder = sys.argv[2] #outfolder
wordn=int(sys.argv[3])
sample=int(sys.argv[4])
threads=int(sys.argv[5])

if __name__ == '__main__':

	#executor = concurrent.futures.ProcessPoolExecutor(threads)
	#futures = [executor.submit(shan, i,wordn) for i in filelist]
	#concurrent.futures.wait(futures)
					
	
	
	for i in filelist:
		shan(i,wordn)

	#print(data)




		