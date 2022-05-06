cd $rnew
o=1
t=2
for i in {0..29}
do
#	echo "$(cat $i$t.fq | wc -l) / 4" | bc
#	cat $i.join.sset | tr -d '\n' | wc -c
	echo "$(cat un.$i.fa |  wc -l) / 2" | bc
done
