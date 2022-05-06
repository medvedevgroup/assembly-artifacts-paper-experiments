o=1
t=2

cd $rnew
for i in {0..29}
do

cat <(fa_to_sset.sh $i.al.fa) <(fa_to_sset.sh $i.un.fa | rev | tr ACGT TGCA) > $i.join.sset
sset_to_fa.sh $i.join.sset > $i.join.fa

bowtie2 -f -x $gloc/refindex.$i  -U $i.join.fa --un un.$i.fa -S $i.join.sam --norc --score-min 'C,0,-1' -p 8

done
