if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi
ref=$1
bam=$2
vcfs=$3
outpref=$4
mpileup=$outpref.mpileup
allelefreqcalls=$outpref.samtools.vcf
filelist=$outpref.filelist.txt

# Run samtools mpileup
if [ ! -r $mpileup ]
then
  samtools mpileup --reference $ref $bam -o $mpileup
fi

# Run samtools-based variant calling
javac $BINDIR/src/*.java
java -cp $BINDIR/src CallVariants pileup_file=$mpileup out_file=$allelefreqcalls

# Print out vcf filenames with absolute paths to filelist
vcfarray=$(echo $vcfs | tr "," "\n")

# Output vcf filenames to a list
if [ -r $filelist ]
then
  rm $filelist
fi
touch $filelist
for vcf in $vcfarray
do
    readlink -f $vcf >> $filelist
done
readlink -f $allelefreqcalls >> $filelist

# Run merging
java -cp $BINDIR/src MergeVariants illumina_bam=None file_list=$filelist out_file=$outpref.all_callers.combined.vcf
# Print possible false positives
#cat $outpref.merged.vcf | grep 'SUPP_VEC=10;' > $outpref.check.txt
# Print possible false negatives
#cat $outpref.merged.vcf | grep 'SUPP_VEC=01;' >> $outpref.check.txt

# Print inconsistent variants when there are 3 samples
#cat $outpref.merged.vcf | grep 'SUPP_VEC=110;' >> $outpref.check.txt
#cat $outpref.merged.vcf | grep 'SUPP_VEC=101;' >> $outpref.check.txt
#cat $outpref.merged.vcf | grep 'SUPP_VEC=100;' >> $outpref.check.txt
#cat $outpref.merged.vcf | grep 'SUPP_VEC=011;' >> $outpref.check.txt
#cat $outpref.merged.vcf | grep 'SUPP_VEC=010;' >> $outpref.check.txt
#cat $outpref.merged.vcf | grep 'SUPP_VEC=001;' >> $outpref.check.txt
