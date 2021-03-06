if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

tablefn=$1
genefn=$2
genomefn=$3
outpref=$4

consensusvcf=$outpref.consensus.vcf
allvcf=$outpref.allsnps.vcf

consensuscombinedvcf=$outpref.consensus.combined.vcf
allcombinedvcf=$outpref.allsnps.combined.vcf

# Compile code
javac $BINDIR/src/*.java

# Convert post-filtering table to VCFs
java -cp $BINDIR/src TableToVcf table_file=$tablefn consensus_file=$consensusvcf all_file=$allvcf

# Combine adjacent SNPs in VCFs
java -cp $BINDIR/src CombineVariants vcf_file=$consensusvcf gene_file=$genefn genome_file=$genomefn out_file=$consensuscombinedvcf
java -cp $BINDIR/src CombineVariants vcf_file=$allvcf gene_file=$genefn genome_file=$genomefn out_file=$allcombinedvcf
