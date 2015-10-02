# Tip and tricks for BAM files - Summary
- [Usefull tools](https://github.com/IARC-bioinfo/BAM-tricks#usefull-tools)
- [Tip and tricks](https://github.com/IARC-bioinfo/BAM-tricks#tip-and-tricks)


## Usefull tools
### Samtools organisation and repositories
- File format [specification](http://samtools.github.io/hts-specs/)
- Samtools [github page](https://github.com/samtools/samtools)
- Samtools [webpage](http://www.htslib.org)
- Samtools [man page](http://www.htslib.org/doc/samtools.html)

Compilation (from [here](http://samtools.github.io/bcftools/)):
```bash
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
git clone --branch=develop git://github.com/samtools/samtools.git
cd bcftools; make
cd ../samtools; make
cd ../htslib; make
```

### Other tools
- [biobambam2](https://github.com/gt1/biobambam2) github
- [Picard](http://broadinstitute.github.io/picard/) with its [github page](https://github.com/broadinstitute/picard)
- R package [rbamtools](https://cran.r-project.org/web/packages/rbamtools/index.html) 
- bedtools [documentation](http://bedtools.readthedocs.org) and [github page](https://github.com/arq5x/bedtools2)

## Tip and tricks
### Check if BAM file is sorted
```bash
samtools view -H test.bam | grep @HD
```
Will show `SO:coordinate` if sorted.
### Sort a BAM file
```bash
samtools sort -o test_sorted.bam -T tmp test.bam
```

### Extract sample name
Only consider the first read group
```bash
samtools view -H test.bam | grep @RG | head -1 | sed "s/.*SM:\([^\t]*\).*/\1/"
```

### Change sample name 
For all read groups in a BAM file
```bash
samtools view -H test.bam  | sed "s/SM:[^\t]*/SM:TEST_SAMPLE_NAME/g" | samtools reheader - test.bam > test_SM.bam
```

### Simple variant calling using freebayes
```bash
freebayes -f ucsc.hg19.fasta --min-alternate-count 5 --no-complex --min-mapping-quality 20 --min-base-quality 20 \
  --min-coverage 20 test.bam | vcffilter -f "QUAL > 20" |  vcfbreakmulti | vt normalize - -q -r ucsc.hg19.fasta > test.vcf 
```

### Add a new ID-specified read group in the bam header 
```bash
RGline=$(samtools view -H test.bam | grep '@RG' | head -1 | sed "s/ID:[^\t]*/ID:NEW_ID/g")
line_number=$(samtools view -H test.bam | sed -n '/@RG/=' - | head -1)
samtools view -H test.bam | sed "${line_number}"'i\'"${RGline}" | samtools reheader - test.bam > test_RGadded.bam
unset RGline; unset line_number
```

### Calculate coverage for each position of a bed
```bash
bedtools coverage -d  -a my_bed.bed -b my_bam.bam
```
