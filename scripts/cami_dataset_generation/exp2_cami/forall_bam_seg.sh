
RLOC=~/s/22/cami/CAMISIM/out_spades_final/2022.03.29_03.02.08_sample_0/reads/
BLOC=~/s/22/cami/CAMISIM/out_spades_final/2022.03.29_03.02.08_sample_0/bam/
GLOC=~/s/22/cami/CAMISIM/spades-genomes

PPATH=$(realpath $(dirname "$0"))
echo $PPATH
k=75


mkdir -p $RLOC/segs$k
rm $RLOC/segs$k/seg.all
touch $RLOC/segs$k/seg.all

cd $RLOC/segs$k
for i in {0..29}
do
	echo "doing seg for $i"
	python $PPATH/segment.py $BLOC/$i.awkresult $GLOC/g.all.fa $k "$RLOC/segs$k/seg.$i"
	cat seg.$i >> seg.all

	sset_to_fa.sh seg.$i > seg.$i.fa
	bowtie2-build -p -f seg.$i.fa segindex.$i

	#bowtie2 -f -x $GLOC/refindex.$i  -U seg.$i.fa -S seg_ref.$k.sam --score-min 'C,0,-1' -p 8
	#igvtools sort seg_ref.$k.sam sorted_seg_ref.$k.sam
	#igvtools index sorted_seg_ref.$k.sam
done


sset_to_fa.sh seg.all > seg.all.fa
bowtie2-build -p -f seg.all.fa segindex
