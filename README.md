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

## CheckVariants

This software validates the presence/absence of SNPs based on SAM-format alignments

### Running

```
Usage: java -cp src CheckVariants [args]
  Example: java -cp src CheckVariants sam_file=jhu004.sam vcf_file=jhu004.vcf genome_file=ref.fa

Required args:
  sam_file    (String) - a SAM file with the alignments of the reads
  vcf_file    (String) - a VCF file with the variant calls
  genome_file (String) - a FASTA file with the reference genome

Optional args:
  coverage_threshold  (int)    [20]    - the min coverage needed to possibly flag a region
  genome_max_len      (int)    [31000] - an upper bound on the genome length
  missed_variant_freq (float)  [0.6]   - call a possible missed variant if ref allele frequency < this value
  fp_freq             (float)  [0.4]   - call a false positive if ref allele frequency > this value
```
