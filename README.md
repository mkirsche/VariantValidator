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

This software merges VCFs from multiple sources.

### Running

```
Usage: java -cp src MergeVariants [args]
  Example: java -cp src MergeVariants filelist=vcflist.txt out_file=merged.vcf

Required args:
  file_list   (String) - a txt file containing absolute paths to VCF files, one on each line
  out_file    (String) - file to write merged variants to
 ```

## CombineVariants

This software takes a merged VCF and combines variants in adjacent positions.

### Running

```
Usage: java -cp src CombineVariants [args]
  Example: java -cp src CombineVariants vcf_file=merged.vcf out_file=merged_combined.vcf

Required args:
  vcf_file    (String) - vcf file containing the variants after merging across samples
  out_file    (String) - file to output variants after combining adjact positions
 ```

## TableToVcf

This software takes a post-filtered TSV and creates 2 VCF files - one with all variants and one with consensus variants

### Running

```
Usage: java -cp src TableToVcf [args]
  Example: java -cp src TableToVcf table_file=merged.vcf consensus_file=consensus.vcf all_file=all.vcf

Required args:
  table_file     (String) - vcf file containing the variants after merging across samples
  consensus_file (String) - file to output consensus variants to
  all_file       (String) - file to output all variants to
 ```


## Recommended Pipeline

Inputs are `sample.fa`, `sample.vcf` from some variant caller, and `sample.bam`.

The variant calling/merging pipeline has been implemented in run.sh which takes the following parameters:

`./run.sh <reference> <bam alignments> <vcf variant calls, different files separated by commas> <output prefix>`


