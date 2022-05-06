#Input: [1] [2] [3] {4} {5}   
#echo "$(realpath $0)"
CODE_DIR=$(dirname "$(realpath $0)")

K=$1    
R=$2 #reference sset (sset to fa numbered),    
O=$3 #output directory 
U=$4 #directed unitig sset (sset to fa numbered)   [make sure k value is right]    
S=$5 #read/segment file in sset format (expect that all strings are substring of reference)


if [[ -z "$K" ]] && [[ -z "$R" ]] && [[ -z "$O" ]] 
	then echo "Error: Missing argument";
	exit 1
fi

RR=$(realpath $R)
OR=$(realpath $O)

mkdir -p $OR
cd $OR

REF_FA="ref.fa"
sset_to_fa.sh $RR > $REF_FA
if [[ -z "$U" ]] 
	U="unid$K.sset"
	then $CODE_DIR/directed_unitigs.sh $REF_FA $K
fi

U_FA="unid$K.fa"
UR=$(realpath $U)
sset_to_fa.sh $UR > $U_FA



U_FA="unid$K.fa"

touch unaligned 

#build index in test.ebwt
~/w/bowtie/bowtie-build -p -f $REF_FA test
NUM_STRING=$(cat $RR | wc -l)
~/w/bowtie/bowtie -s -a -k $NUM_STRING -v 0 test -f $U_FA --un unaligned --norc > test.sam

NUM_UNITIG=$(cat $UR | wc -l) 
NUM_TRUE_INVALID_UNITIG=$($(cat unaligned | wc -l) / 2)

echo "n_string n_unitig n_true_invalid"
echo "$NUM_STRING $NUM_UNITIG $NUM_TRUE_INVALID_UNITIG"
