# VariantValidator

This software validates the presence/absence of SNPs based on SAM-format alignments

## Compilation

```
javac src/*.java
```

## Running

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
