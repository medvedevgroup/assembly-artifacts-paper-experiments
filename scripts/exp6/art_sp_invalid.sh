if [[ -z $1 || -z $2 || -z $3 || -z $4 ]] ;
then
        echo "required args: <data-parent> <k> <cov> <read-len>" 
        exit
fi



BASEDIR=$(dirname "$(realpath $0)")

D="$(realpath $1)"
K=$2
C=$3


T=20
L=$4

KM1=$(($K-1))

mkdir -p $D/$K/E/$C/$L


WDIR=$D/$K/E/$C/$L
cd $WDIR
mkdir -p tmp_spades


##get simulated reads ~/w/artsimulator...
EF_READ="$D/art_reads/ef_$L.$C.fa"
AC_READ="$D/art_reads/ref_$L.$C.sam.fq"


if [[ 1 -eq 0 ]] ; then

### GATB_UNITIG ####
GATB_UNI=ef_$L.$C.unitigs.fa


cat unitig_gatb_to_ref_sort.sam | awk '$3=="contig_0" {print $10, length($10), $4}' > gatb_uni_align.txt 



bowtie2 -f -x $D/refindex_bt2  -U $WDIR/$GATB_UNI -S unitig_gatb_to_ref.sam --score-min 'C,0,-1' -p $T
bowtie2 -f -x $D/refindex_bt2  -U $WDIR/$GATB_UNI -S unitig_gatb_to_ref_indel.sam -p $T

igvtools sort unitig_gatb_to_ref.sam unitig_gatb_to_ref_sort.sam
igvtools index unitig_gatb_to_ref_sort.sam

igvtools sort unitig_gatb_to_ref_indel.sam unitig_gatb_to_ref_indel_sort.sam
igvtools index unitig_gatb_to_ref_indel_sort.sam


##### ABYSSS ########
#cp $D/$K/r/palin_align.txt $WDIR/

~/w/abyss/bin/ABYSS $EF_READ -k $K -o $WDIR/contig_abyss.fa
#cp $BASEDIR/contig_abyss.fa $WDIR/contig_abyss.fa

bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_abyss.fa --un invalid_con_abyss -S contig_abyss_to_ref.sam --score-min 'C,0,-1' -p $T
bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_abyss.fa --un invalid_con_abyss_indel -S contig_abyss_to_ref_indel.sam -p $T

igvtools sort contig_abyss_to_ref.sam contig_abyss_to_ref_sort.sam
igvtools index contig_abyss_to_ref_sort.sam

igvtools sort contig_abyss_to_ref_indel.sam contig_abyss_to_ref_indel_sort.sam
igvtools index contig_abyss_to_ref_indel_sort.sam


cat contig_abyss_to_ref_sort.sam | awk '$3=="contig_0" {print $10, length($10), $4}' > abyss_align.txt 





###### GATB ###########
mkdir -p tmp_gatb
cd tmp_gatb
~/w/gatb-minia-pipeline/gatb -s $EF_READ  --kmer-sizes $K 
#~/w/gatb-minia-pipeline/gatb -s $EF_READ --no-error-correction --kmer-sizes $K --abundance-mins 1
#~/w/gatb-minia-pipeline/gatb -s $EF_READ
multifa_to_fa.sh assembly.fasta > $WDIR/contig_gatb.fa
rm *.glue*
cd $WDIR

bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_gatb.fa --un invalid_con_gatb -S contig_gatb_to_ref.sam --score-min 'C,0,-1' -p $T
bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_gatb.fa --un invalid_con_gatb_indel -S contig_gatb_to_ref_indel.sam -p $T

igvtools sort contig_gatb_to_ref.sam contig_gatb_to_ref_sort.sam
igvtools index contig_gatb_to_ref_sort.sam

igvtools sort contig_gatb_to_ref_indel.sam contig_gatb_to_ref_indel_sort.sam
igvtools index contig_gatb_to_ref_indel_sort.sam

cat contig_gatb_to_ref_sort.sam | awk '$3=="contig_0" {print $10, length($10), $4}' > gatb_align.txt 



########### SPADES #############
##does error correction
#python ~/w/SPAdes-3.15.3-Linux/bin/spades.py -s $AC_READ -k $KM1 -o tmp_spades -t $T

python ~/w/SPAdes-3.15.3-Linux/bin/spades.py -s $EF_READ -k $KM1 -o tmp_spades -t $T --only-assembler 
multifa_to_fa.sh $WDIR/tmp_spades/contigs.fasta > $WDIR/contig_spades.fa

bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_spades.fa --un invalid_con_spades -S contig_spades_to_ref.sam --score-min 'C,0,-1' -p $T
bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_spades.fa --un invalid_con_spades_indel -S contig_spades_to_ref_indel.sam -p $T 

igvtools sort contig_spades_to_ref.sam contig_spades_to_ref_sort.sam
igvtools index contig_spades_to_ref_sort.sam

igvtools sort contig_spades_to_ref_indel.sam contig_spades_to_ref_indel_sort.sam
igvtools index contig_spades_to_ref_indel_sort.sam

cat contig_spades_to_ref_sort.sam | awk '$3=="contig_0" {print $10, length($10), $4}' > spades_align.txt 


fi

####### MEGA #################
rm -rf tmp_mega
~/w/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -r $EF_READ --k-list $KM1  --min-count 1 -o tmp_mega --keep-tmp-files
multifa_to_fa.sh tmp_mega/final.contigs.fa > $WDIR/contig_mega.fa

#GET ONLY INVALIDS MEGAHIT
bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_mega.fa --un invalid_con_mega -S contig_mega_to_ref.sam --score-min 'C,0,-1' -p $T
bowtie2 -f -x $D/refindex_bt2  -U $WDIR/contig_mega.fa --un invalid_con_mega_indel -S contig_mega_to_ref_indel.sam -p $T


igvtools sort contig_mega_to_ref.sam contig_mega_to_ref_sort.sam
igvtools index contig_mega_to_ref_sort.sam

igvtools sort contig_mega_to_ref_indel.sam contig_mega_to_ref_indel_sort.sam
igvtools index contig_mega_to_ref_indel_sort.sam

cat contig_mega_to_ref_sort.sam | awk '$3=="contig_0" {print $10, length($10), $4}' > mega_align.txt 
exit 1



#~/w/SPAdes-3.15.3-Linux/bin/spades-gbuilder $AC_READ $WDIR/spades_unitig_multi.fa -k $KM1 -t $T
#~/w/SPAdes-3.15.3-Linux/bin/spades-gbuilder $EF_READ $WDIR/spades_unitig_ef_multi.fa -k $KM1 -t $T
#python ~/w/SPAdes-3.15.3-Linux/bin/spades.py -s $AC_READ -k $KM1 -o tmp_spades --only-error-correction


#multifa_to_fa.sh $WDIR/spades_unitig_multi.fa > $WDIR/spades_unitig.fa
#multifa_to_fa.sh $WDIR/spades_unitig_ef_multi.fa > $WDIR/spades_unitig_ef.fa


####~/w/bowtie/bowtie-build -p -f $D/ref.fa refindex ### index was built like this
#~/w/bowtie/bowtie -s -a -k 2 -v 0 $D/refindex -f $WDIR/spades_unitig.fa --un invalid_spades -p $T > spades_invalid_unitig.sam
#~/w/bowtie/bowtie -s -a -k 2 -v 0 $D/refindex -f $WDIR/spades_unitig_ef.fa --un invalid_spades_ef -p $T > spades_invalid_unitig_ef.sam

#bowtie -s -a -k 2 -v 0 $D/refindex_bt2 -f $WDIR/spades_contig.fa --un invalid_contigs_spades   > spades_invalid_contig.sam
 


#cat contig_abyss_to_ref_sort.sam | awk '$3=="NC_000024.10" {print $10, length($10), $4}' > abyss_align.txt 
#cat contig_mega_to_ref_sort.sam | awk '$3=="NC_000024.10" {print $10, length($10), $4}' > mega_align.txt 
#cat contig_spades_to_ref_sort.sam | awk '$3=="NC_000024.10" {print $10, length($10), $4}' > spades_align.txt 
#cat contig_gatb_to_ref_sort.sam | awk '$3=="NC_000024.10" {print $10, length($10), $4}' > gatb_align.txt 
