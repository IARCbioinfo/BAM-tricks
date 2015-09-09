# Tip and tricks for BAM files

## Usefull tools
### Samtools organisation and repositories
- File format specification: http://samtools.github.io/hts-specs/
- Samtools github page https://github.com/samtools/samtools
- Samtools webpage http://www.htslib.org
- Samtools man page: http://www.htslib.org/doc/samtools.html

### Other tools
- [biobambam2](https://github.com/gt1/biobambam2) github
- [Picard](http://broadinstitute.github.io/picard/) with its [github page](https://github.com/broadinstitute/picard)
- R package [rbamtools](https://cran.r-project.org/web/packages/rbamtools/index.html) 

## Tip and tricks
### Check if BAM file is sorted (`SO:coordinate`)
```bash
samtools view -H test.bam | grep @HD
```

### Sort a BAM file
```bash
samtools sort -o test_sorted.bam -T tmp test.bam
```

### Extract sample name from a BAM file using samtools (only consider the first read group)
```bash
samtools view -H test.bam | grep @RG | head -1 | sed "s/.*SM:\([^\t]*\).*/\1/"
```

### Change sample name for all read groups in a BAM file
```bash
samtools view -H test.bam  | sed "s/SM:[^\t]*/SM:TEST_SAMPLE_NAME/g" | samtools reheader - test.bam > test_SM.bam
```

### Do a simple variant calling using freebayes
```bash
freebayes -f ucsc.hg19.fasta --min-alternate-count 5 --no-complex --min-mapping-quality 20 --min-base-quality 20 --min-coverage 20 test.bam | vcffilter -f "QUAL > 20" |  vcfbreakmulti | vt normalize - -q -r ucsc.hg19.fasta > test.vcf 
```
