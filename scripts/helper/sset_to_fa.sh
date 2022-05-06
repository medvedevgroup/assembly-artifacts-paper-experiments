#BAZE=$(basename $1 .sset)
cat $1 | awk '{print ">contig_" NR-1 "\n" $0}'
