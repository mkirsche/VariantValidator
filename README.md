# VariantValidator

Software for calling small variants (or validating existing calls) in a small genome

## Compilation (all software)

```
javac src/*.java

```

## CallVariants

This software uses the output of samtools mpileup to call variants based on simple allele frequency thresholds.

### Running

```
Usage: java -cp src CallVariants [args]
  Example: java -cp src CallVariants pileup_file=jhu004.mpileup out_file=calls.vcf

Required args:
  pileup_file (String) - the output of samtools mpileup with the read alignments and reference
  out_file    (String) - the output file to call SNPs

Optional args:
  coverage_threshold  (int)    [20]    - the min coverage needed to possibly flag a region
  genome_max_len      (int)    [31000] - an upper bound on the genome length
  alt_threshold       (float)  [0.15]   - call a variant if any alt allele frequency > this value
  ref_threshold       (float)  [0.60]   - call an N even if no alt allele frequency is high enough if ref allele frequency < this value)
```

## MergeVariants

This software merges VCFs from multiple sources

### Running

```
Usage: java -cp src MergeVariants [args]
  Example: java -cp src MergeVariants filelist=vcflist.txt out_file=merged.vcf

Required args:
  file_list   (String) - a txt file containing absolute paths to VCF files, one on each line
  out_file    (String) - file to write merged variants to
 ```

## Recommended Pipeline

Inputs are `sample.fa`, `sample.vcf` from some variant caller, and `sample.bam`.

```
samtools mpileup --reference sample.fa sample.bam -o sample.mpileup
javac src/*.java
java -cp src CallVariants pileup_file=sample.mpileup out_file=allele_freq_calls.vcf
# Print out vcf filenames with absolute paths to filelist.txt
readlink -f sample.vcf > filelist.txt
readlink -f allele_freq_calls.vcf >> filelist.txt
java -cp src MergeVariants file_list=filelist.txt out_file=merged.vcf
# Print possible false positives
cat merged.vcf | grep 'SUPP_VEC=10'
# Print possible false negatives
cat merged.vcf | grep 'SUPP_VEC=01'
```
