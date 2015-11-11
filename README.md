# Tip and tricks for BAM files
- [Usefull tools](https://github.com/IARC-bioinfo/BAM-tricks#usefull-tools)
 - [Samtools organisation and repositories](https://github.com/IARC-bioinfo/BAM-tricks#samtools-organisation-and-repositories) 
 - [Other tools](https://github.com/IARC-bioinfo/BAM-tricks#other-tools)
- [Tip and tricks](https://github.com/IARC-bioinfo/BAM-tricks#tip-and-tricks)
 - [Check if BAM is sorted](https://github.com/IARC-bioinfo/BAM-tricks#check-if-bam-file-is-sorted) 
 - [Sort a BAM](https://github.com/IARC-bioinfo/BAM-tricks#sort-a-bam-file)
 - [Extract sample name](https://github.com/IARC-bioinfo/BAM-tricks#extract-sample-name)
 - [Change sample name](https://github.com/IARC-bioinfo/BAM-tricks#change-sample-name)
 - [Simple variant calling with freebayes](https://github.com/IARC-bioinfo/BAM-tricks#simple-variant-calling-using-freebayes)
 - [Add new read group in header](https://github.com/IARC-bioinfo/BAM-tricks#add-a-new-id-specified-read-group-in-the-bam-header)
 - [Calculate bed-positions coverage](https://github.com/IARC-bioinfo/BAM-tricks#calculate-coverage-for-each-position-of-a-bed)


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
bedtools coverage -d -a my_bed.bed -b my_bam.bam
```

### Find where reads map in a BAM file

This will only report non-zero coverage and create contiguous regions with similar coverage (see [`bedtools` manual](http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html)). Here for example we only keep regions with at least 10 reads:
```bash
$ bedtools genomecov -bg -ibam test.bam | awk '$4>=10'
chr1	154566162	154566163	14
chr1	154566163	154566164	15
chr1	154566164	154566167	18
chr1	154566167	154566171	19
```

It can also be useful to group these regions and to report the total number of reads in each region. The first step is to merge contiguous regions identified by the previous command with [`bedtools merge`](http://bedtools.readthedocs.org/en/latest/content/tools/merge.html). Here we also merge regions less than 100bp from each other:
```bash
$ bedtools genomecov -bg -ibam test.bam | awk '$4>=10' | bedtools merge -d 100 -i stdin
chr1	154566162	154566294
chr2	45189642	45189787
```

Finally we can now ask `bedtools` to count the number of reads in each of these regions using [`coverage`](http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html), and the full command is:
```bash
$ bedtools genomecov -bg -ibam test.bam | awk '$4>=10' | bedtools merge -d 100 -i stdin | bedtools coverage -a stdin -b test.bam
chr1	154566162	154566294	21	132	132	1.0000000
```
Here for example there are 21 reads in region chr1:154566162-154566294. These 21 reads cover 132bp and the region itself is 132bp long, so 100% of the region is covered.
