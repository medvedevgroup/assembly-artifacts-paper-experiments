k=75
#for k in {21,31,41,51,61,71}
for c in 1
do
        cd $rnew

       jellyfish count -m $k -s 100M -t 10 all.fa -o mer_counts.jf
        python $rloc/directed/newseg.py mer_counts.jf $gloc/ref.all.fa $k seg.all.sset
#        sset_to_fa.sh  seg_test.sset > seg_test.fa

        bowtie2-build -p -f seg.all.fa segindex

#       sset_to_fa.sh unid$k.sset > unid$k.fa
#        bowtie2 -f -x segtest  -U unid$k.fa -S unid$k.test.sam --un unsafe.test.fa --score-min 'C,0,-1' -p 8 --norc
#       bowtie2 -f -x segnew  -U unid$k.fa -S unid$k.sam --un unsafe.fa --score-min 'C,0,-1' -p 8 --norc
#       bowtie2 -f -x ../../../refindex  -U unid$k.fa -S unid$k.mis.sam --un misassembled.fa --score-min 'C,0,-1' -p 8 --norc

#       echo "$k $c"
        echo "$(cat unsafe.test.fa | wc -l) / 2" | bc
#       echo "$(cat misassembled.fa | wc -l) / 2" | bc
#        cd ../../../

done
