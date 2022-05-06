BASEDIR=$(dirname "$(realpath $0)")
#echo "$BASEDIR"
D=$(realpath $1)
K=$2
S=$3
mkdir -p $D/$K/d/$S
cd $D/$K/d/$S

#s gap len
#T gap freq
python $BASEDIR/simu_read_rand.py $D/ref.sset $D/$K/d/$S/ref_read.sset $S $K

sset_to_fa.sh $D/$K/d/$S/ref_read.sset > $D/$K/d/$S/ref_read.fa
$BASEDIR/dd_unitigs.sh $D/$K/d/$S/ref_read.fa $K


REF_FA=$D/ref.fa

#GET ONLY INVALIDS
cd $D/$K/d/$S/
sset_to_fa.sh $D/$K/d/$S/unid$K.sset > unid$K.fa

#build index in refindex.ebwt
#~/w/bowtie/bowtie-build -p -f $REF_FA refindex
NUM_STRING=$(cat $D/ref.sset | wc -l)
#~/w/bowtie/bowtie -s -a -k $NUM_STRING -v 0 refindex -f unid$K.fa --un invalids_true > ref_uni.sam
~/w/bowtie/bowtie -a -k $NUM_STRING -v 0 refindex -f unid$K.fa --un invalids_true > ref_uni.sam
python $BASEDIR/keep_canonical.py invalids_true invalid_bcalm_canon


if false; then
sset_to_fa.sh segfile.sset > segfile.fa
~/w/bowtie/bowtie-build -p -f segfile.fa segindex
NUM_SEG=$(cat segfile.sset | wc -l)
~/w/bowtie/bowtie -s -a -k $NUM_SEG -v 0 segindex -f unid$K.fa --un invalids > seg_uni.sam



echo "reporting: "
N_TINVALID=$($(cat invalids_true | wc -l) / 2)
N_INVALID=$($(cat invalids | wc -l) / 2)
N_UNI=$(wc -l unid$K.sset | cut -f1 -d" " )
N_KMER=$(cat n_mer$K)
N_READ=$($(cat ref_read.sset | wc -l) / 2)
N_SEG=$($(cat segfile.sset | wc -l) / 2)

echo "k coverage n_unitig n_invalid n_true_invalid n_kmer n_read n_seg n_true_invalid_by_spades"
echo "$K $S $N_UNI $N_INVALID $N_TINVALID $N_KMER $N_READ $N_SEG" > report.txt
cat report.txt
fi


#GET TRULY INVALID UNITIGS





### ALL 
#$BASEDIR/pipe_pseudo_read.sh $D $K $S
