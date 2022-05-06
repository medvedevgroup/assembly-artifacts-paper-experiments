
import re
import sys
import math
import os

fileGenome = sys.argv[1]  # 5 columns
fileSeg = sys.argv[2] # multifa to fa
fileOut=(sys.argv[3])


genome_ids=[]
genome_id_to_seq={}

fixed_gid_to_loc={}


reads=[]

total_bp_all=0
num_reads_all=0

with open(fileGenome) as f: #read reference genome fasta, g.10.fa
    gcount=0
    id=""
    seq=""
    while True:
        gcount += 1
        g = f.readline()
        
        # if line is empty # end of file is reached
        if not g:
            break
        
        if g[0]=='>':
            id = g.split(" ")[0][1:].rstrip().lstrip().lower()
            genome_ids.append(id)
        else:
            seq=g.rstrip()
            genome_id_to_seq[id]=seq        

read_glocs=[]
read_gids=[]


wf = open(fileOut, 'w')

ref_seq=genome_id_to_seq["contig_0"]
with open(fileSeg) as f: #read reference genome fasta, g.10.fa
    gcount=0
    id=""
    seq=""
    while True:
        gcount += 1
        g = f.readline()
        
        # if line is empty # end of file is reached
        if not g:
            break
        
        if g[0]!='>':
            seq=g.rstrip().lstrip()
            wf.write(str(ref_seq.find(seq))+"\n")

wf.close()
print("v3_num_seg=",numseg)


