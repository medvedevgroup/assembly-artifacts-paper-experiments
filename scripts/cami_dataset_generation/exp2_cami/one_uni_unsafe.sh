
k=75
cd  $rnew/d_all_$k
#bowtie2 -a -f -x $sloc/segs$k/segindex  -U unid$k.fa -S uni_seg.all.sam --un unsafe.fa --score-min 'C,0,-1' -p 8
sset_to_fa.sh unid$k.sset > unid$k.fa
bowtie2 -a -f -x $gloc/refindex.all  -U unid$k.fa -S uni_ref.all.sam --un all_mis.fa --score-min 'C,0,-1' -p 8

