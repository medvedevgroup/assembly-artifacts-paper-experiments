
# gap_len = 5
# gap_per = 1/1000
import random
import sys
import math

random.seed(1)

genomessetFile = sys.argv[1]
outFile = sys.argv[2]
gap_len = int(sys.argv[3])
gap_freq = int(sys.argv[4])
gap_per=1.0/gap_freq

wf = open(outFile,"w")

with open(genomessetFile) as f:
    g=f.read()
    glist = list(g)
    total_bp = len(g)
    gap_in = int(1/gap_per)

    ranlist_ind=random.sample(range(0, total_bp), math.floor(total_bp*gap_per*1.0))
    #print(ranlist_ind)

    for j in ranlist_ind:
        for m in range(gap_len):
            if(j+m<total_bp):
                glist[j+m] = 'N'


    # j = gap_in
    # while True:
    #     for m in range(gap_len):
    #         if(j+m<total_bp):
    #             glist[j+m] = 'N'
    #     j = j + gap_len +  gap_in
    #     if j >= total_bp:
    #         break
    g="".join(glist)

    wf.write(g)
    wf.close()
    print("total_bp=",total_bp)
    print("gap_len=",gap_len)
    print("gap_per=",gap_per)
    print("N_count = ", g.count('N') )
