import sys
import random
import math

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
if(len(sys.argv)<2):
    print("required args: <fasta/sset ref>")
    sys.exit(1)

filepath=sys.argv[1]
ref=""
with open(filepath) as fp:
   line = fp.readline().rstrip()
   while line:
       if(line[0]!='>'):
           ref = line
       line = fp.readline().rstrip()

#sys.exit(2)
# refl=["AC"]*100000
# ref=''.join(refl)


total_bp=len(ref)
num_palin=math.floor(total_bp/10000)
palin_len_half=500
palin_len=palin_len_half*2


ranlist_ind=random.sample(range(0, total_bp-palin_len), num_palin)
    #print(ranlist_ind)


ref_chars = list(ref)

for j in ranlist_ind:
    ii=0
    for i in range(palin_len-1, palin_len_half-1, -1):
        ref_chars[j+i] = complement[ref_chars[j+ii]]
        ii+=1

newref=''.join(ref_chars)

print(newref)
