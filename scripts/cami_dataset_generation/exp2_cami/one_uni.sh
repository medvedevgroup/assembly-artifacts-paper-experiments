k=55

cd $rnew
#cat *.join.sset > all.sset
#sset_to_fa.sh all.sset > all.fa

mkdir -p d_all_$k
cd d_all_$k
dsk  -file ../all.fa  -kmer-size $k -nb-cores 8 -out all_$k.h5 -abundance-min 1
dsk2ascii -file all_$k.h5 -out all_$k.ascii -nb-cores 8 
sset_to_fa.sh <(cat all_$k.ascii | cut -f 1 -d" ") > all_$k.fa 
~/s/21/unitig-validator/directed_unitigs.sh all_$k.fa $k 
