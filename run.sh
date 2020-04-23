BINDIR=`pwd`
ref=$1
bam=$2
vcf=$3
outpref=$4
mpileup=$outpref.mpileup
allelefreqcalls=$outpref.allele_freq_calls.vcf
filelist=$outpref.filelist.txt
if [ ! -r $mpileup ]
then
  samtools mpileup --reference $ref $bam -o $mpileup
fi
javac $BINDIR/src/*.java
java -cp $BINDIR/src CallVariants pileup_file=$mpileup out_file=$allelefreqcalls
# Print out vcf filenames with absolute paths to filelist
readlink -f $vcf > $filelist
readlink -f $allelefreqcalls >> $filelist
java -cp $BINDIR/src MergeVariants file_list=$filelist out_file=$outpref.merged.vcf
# Print possible false positives
cat $outpref.merged.vcf | grep 'SUPP_VEC=10' > $outpref.check.txt
# Print possible false negatives
cat $outpref.merged.vcf | grep 'SUPP_VEC=01' >> $outpref.check.txt
