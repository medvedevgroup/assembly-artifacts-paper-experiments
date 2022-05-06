o=1
t=2

cd $rnew
for i in {0..29}
do


cat $i$o.fq $i$t.fq > $i.fq
fq_to_fa.sh $i.fq > $i.fa 
bowtie2 -f -x $gloc/refindex.$i  -U $i.fa -S $i.sam --un $i.un.fa --al $i.al.fa --score-min 'C,0,-1' -p 8 --norc

done
