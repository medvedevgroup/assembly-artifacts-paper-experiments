
import re
import sys
import math
import os

fileMer = sys.argv[1]  # 5 columns
fileGenome = sys.argv[2] # multifa to fa
k=int(sys.argv[3])
fileSeg = sys.argv[4] 


# fileSamCut = "/Users/aur1111/projects/server-sync/blightdir/20.awkresult"
# fileGenome = "/Users/aur1111/projects/server-sync/g.20.fa"
# k=55
# fileSeg = "seg"


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

# for i in genome_id_to_seq:
#     print(i)
numseg=0
wf = open(fileSeg,"w")
wf2 = open("temp","w")
for gid in genome_ids:
    seq = genome_id_to_seq[gid]
    for i in range(len(seq) - k + 1):
        kmer=seq[i:i + k]
        wf2.write(kmer+"\n")
wf2.close()



os.system("cat temp | jellyfish query -i {0} | cut -f2 -d \" \" > outkmer".format(fileMer))

file1 = open("outkmer", 'r')
Lines = file1.readlines()         

linecount=0
for gid in genome_ids:
    seq = genome_id_to_seq[gid]
    segstart=0
    segend=-1
    seg=""
    for i in range(len(seq) - k + 1):
        kmer=seq[i:i + k]
        #check if kmer in jellyfish
        #os.system("cat {0} | cut -f2 -d\" \" | jellyfish query -i mer_counts.jf  | cut -f2 -d \" \" > {1}".format(vertex_f, copycount_f))
        count=int(Lines[linecount].strip())
        linecount+=1
        if(count==0):
            if(seg!=""):
                wf.write(seg+"\n")
#                print(segstart, segend)		
                numseg+=1
                seg=""
        else:
            if seg=="":
                segstart=i
                seg = seg+kmer
                segend = i+k-1
            else:
                seg = seg+kmer[-1]
                segend = segend + 1
    if seg!="":
        wf.write(seg+"\n")
#        print(segstart, segend)  
        numseg+=1       
  

wf.close()
print("v3_num_seg=",numseg)


