
import re
import sys
import math

fileSamCut = sys.argv[1]  # 5 columns
fileGenome = sys.argv[2] # multifa to fa
k=int(sys.argv[3])
fileSeg = sys.argv[4] 

# fileSamCut = "/Users/aur1111/projects/server-sync/blightdir/20.awkresult"
# fileGenome = "/Users/aur1111/projects/server-sync/g.20.fa"
# k=55
# fileSeg = "seg"

read_len = k

genome_ids=[]
genome_id_to_seq={}

fixed_gid_to_loc={}

# samcutfile = "$i.samcut"
# outFile = "segfile.sset"
# readlen = 250

wf_readloc = open("readloc","w")

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


fff=0
lll=1
def first_key(n):
    return n[fff]  

def last_key(n):
    return n[lll]  

read_glocs_dup={}

with open(fileSamCut) as f:
    gcount=0
    while True:
        gcount += 1

        # Get next line from file
        line = f.readline()

        # if line is empty
        # end of file is reached
        if not line:
            break

        line=line.rstrip()
        
        read_gid = line.split()[0].rstrip().lstrip().lower()
        read_gloc=int(line.split()[1])-1
        read_string=line.split()[2]
        read_len=int(line.split()[3])
#        if(len(read_string)!=read_len):
#            print("exc read len", len(read_string))
        read_glocs_dup[(read_gloc, read_gid)]=""
        # if(read_gid=="anla01000000.2"):
        #     print(str(read_gloc)+" "+str(read_gid)+"\n")
            
        if read_gid not in fixed_gid_to_loc:
            fixed_gid_to_loc[read_gid] = set()
        fixed_gid_to_loc[read_gid].add(read_gloc)


        #read_gids.append(read_gid)
        # if read_gid in read_gid_to_str: 
        #     read_gid_to_str[read_gid].append(line.split()[1])
        # else:
        #     read_gid_to_str[read_gid] = []

wf = open(fileSeg,"w")
 #file id, loc(0 based), string, read len
        # read all indices from file 
   
   
read_glocs=[]     
for tup in read_glocs_dup:
    read_glocs.append(tup)
    
a=sorted(read_glocs, key = first_key)
#b=sorted(a, key = last_key)
read_glocs=a


curr_id=""
gid_glocs={}
for tup in read_glocs:
    read_gloc = tup[0]
    read_gid = tup[1]
    if (curr_id!=read_gid):
        curr_id=read_gid
        gid_glocs[curr_id]=[]
        gid_glocs[curr_id].append(read_gloc)
        #wf_readloc.write(str(read_gloc)+" "+str(read_gid)+"\n")
    else:
        gid_glocs[curr_id].append(read_gloc)
        #wf_readloc.write(str(read_gloc)+" "+str(read_gid)+"\n")

for gid in genome_ids:
    
    if gid not in gid_glocs:
        gid_glocs[gid] = []
    wf_readloc.write(str(gid)+"\n")
    for gloc in gid_glocs[gid]:
        
        wf_readloc.write(str(gloc)+" "+str(id)+"\n")
                
#sys.exit() 
numseg=0
for gid in genome_ids:
    wf_readloc.write(str(gid)+"\n")
    # for gloc in gid_glocs[gid]:
    #     wf_readloc.write(str(gloc)+" "+str(id)+"\n")
    g=genome_id_to_seq[gid]
    
    if gid not in fixed_gid_to_loc:
        continue
    
    read_glocs=sorted(list(fixed_gid_to_loc[gid]))
    # for i in (read_glocs):
    #     wf_readloc.write(str(i))
    # sys.exit()
    
    if(len(fixed_gid_to_loc[gid])!=0):
        segs=[]
        same=False
        segstart=read_glocs[0]
        segend=read_glocs[0]+read_len-1
        for i in range(1,len(read_glocs)):
            
            if(segend - read_glocs[i] + 1 < k-1):
                segs.append(g[segstart:segend+1])
                #wf.write(g[segstart:segend+1]+" "+str(gid)+" "+str(segstart)+" "+str(segend)+" "+str(segend-segstart)+" "+str(len(g))+"\n")
                wf.write(g[segstart:segend+1]+"\n")

                
                #print(segstart, segend,  segend-segstart)
                wf_readloc.write(str(segstart)+" "+str(segend)+" "+str(segend-segstart)+"\n")
                segstart=read_glocs[i]

                same=False
                #segend=ranlist_ind[i]+read_len-1
            else:
                same=True

            segend = read_glocs[i]+read_len-1
        segs.append(g[segstart:segend+1])
        #wf.write(g[segstart:segend+1]+" "+str(gid)+" "+str(segstart)+" "+str(segend)+" "+str(segend-segstart)+" "+str(len(g))+"\n")
        #print(segstart, segend, segend-segstart)
        wf.write(g[segstart:segend+1]+"\n")

        for s in segs:
            #wf.write(s+"\n")
            numseg+=1

        
# for gid in gid_glocs:
#     read_glocs = gid_glocs[gid]
#     g=genome_id_to_seq[gid]
#     if(len(read_glocs)!=0):
#         segs=[]
#         same=False
#         segstart=read_glocs[0]
#         segend=read_glocs[0]+read_len-1
#         for i in range(1,len(read_glocs)):
#             if(segend - read_glocs[i] + 1 < k-1):
#                 segs.append(g[segstart:segend+1])
#                 #print(segstart, segend,  segend-segstart)
#                 wf_readloc.write(str(segstart)+" "+str(segend)+"\n")
#                 segstart=read_glocs[i]
#                 same=False
#                 #segend=ranlist_ind[i]+read_len-1
#             else:
#                 same=True

#             segend = read_glocs[i]+read_len-1
#         segs.append(g[segstart:segend+1])
#         #print(segstart, segend, segend-segstart)

#         for s in segs:
#             wf.write(s+"\n")
#             numseg+=1


wf.close()
print("v3_num_seg=",numseg)

