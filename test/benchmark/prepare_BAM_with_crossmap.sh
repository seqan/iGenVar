# What is CrossMap?
# CrossMap is a program for genome coordinates conversion between different assemblies
# (such as hg18 (NCBI36) <=> hg19 (GRCh37)). It supports commonly used file formats including
# BAM, CRAM, SAM, Wiggle, BigWig, BED, GFF, GTF, MAF VCF, and gVCF.
# Source: http://crossmap.sourceforge.net/, last checked: 06.04.2022.

conda create -n crossmap -c defaults -c bioconda crossmap # Python 2.7.18
# download Chain file:
# download to data
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

echo "$(tput setaf 1)$(tput setab 7)------- CrossMap installed and prepared (5.1/9.6) --------$(tput sgr 0)" 1>&3

# Illumina Mate Pair
## prepare reference file:
## reorder reference:
samtools faidx data/reference/GRCh37/hg19.fa $(cat Repos/iGenVar/test/benchmark/data/order_of_hg19.txt) \
    > data/reference/GRCh37/hg19.reordered.fa
## for GRIDSS:
less data/reference/GRCh37/hg19.reordered.fa | sed -n -e '/>chrM/,$p' | sed -e '/>chr1_gl000191_random/,$d' > data/reference/GRCh37/hg19.reordered.short.fa
less data/reference/GRCh37/hg19.reordered.fa | sed -e '/>chrM/,$d' >> data/reference/GRCh37/hg19.reordered.short.fa

echo "$(tput setaf 1)$(tput setab 7)------- reference files prepared (5.2/9.6) --------$(tput sgr 0)" 1>&3

# Illumina Paired End
samtools sort short_reads/GRCh37/HG002.hs37d5.2x250.bam -o short_reads/GRCh37/HG002.hs37d5.2x250.sorted.bam
samtools sort short_reads/GRCh38/HG002.GRCh38.2x250.bam -o short_reads/GRCh38/HG002.GRCh38.2x250.sorted.bam
samtools index short_reads/GRCh37/HG002.hs37d5.2x250.sorted.bam
samtools index short_reads/GRCh38/HG002.GRCh38.2x250.sorted.bam

# MtSinai PacBio
## Add MD tag for Sniffles
samtools calmd --threads 16 -Q -b long_reads/GRCh37/HG002_PacBio_GRCh37.bam \
    > long_reads/GRCh37/HG002_PacBio_GRCh37.md.bam reference/GRCh37/hs37d5.fa
samtools calmd --threads 16 -Q -b long_reads/GRCh38/HG002_PacBio_GRCh38.bam \
    > long_reads/GRCh38/HG002_PacBio_GRCh38.md.bam reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# PacBio CCS 10kb
## Add MD tag for Sniffles
samtools calmd --threads 16 -Q -b long_reads/GRCh37/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam \
    > long_reads/GRCh37/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.md.bam \
    reference/GRCh37/hs37d5.fa
samtools calmd --threads 16 -Q -b long_reads/GRCh38/HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam \
    > long_reads/GRCh38/HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.md.bam \
    reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
## Add sampled files with different coverages
#   coverage calculation:
#   samtools depth -a --threads 16 DATASET | awk '{sum+=$3} END { print "Average = ",sum/NR}'
cd long_reads/GRCh37/
mkdir -p sampled
### The library was sequenced to approximately 30-fold coverage.
### Samtools sampling is working when sampling 50% or lower:
# 1x coverage
samtools view -s 0.033333333333333333 -b HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam \
    > sampled/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.1x.bam
# 2x coverage
samtools view -s 0.066666666666666666 -b HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam \
    > sampled/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.2x.bam
# 3x coverage
samtools view -s 0.10 -b HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam \
    > sampled/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.3x.bam
# 5x coverage
samtools view -s 0.166666666666666666 -b HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam \
    > sampled/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.5x.bam
# 10x coverage
samtools view -s 0.333333333333333333 -b HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam \
    > sampled/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.10x.bam
cd ../..

# 10X Genomics
## Add MD tag for Sniffles
samtools calmd -b long_reads/GRCh37/NA24385_phased_possorted_bam.bam \
    > long_reads/GRCh37/NA24385_phased_possorted_bam.md.bam \
    reference/GRCh37/hg19.reordered.fa
## Convert from GRCh37 to GRCh38
conda activate crossmap
CrossMap.py bam hg19ToHg38.over.chain.gz long_reads/GRCh37/NA24385_phased_possorted_bam.bam \
    long_reads/GRCh38/NA24385_phased_possorted_bam.Hg38.bam
CrossMap.py bam hg19ToHg38.over.chain.gz long_reads/GRCh37/NA24385_phased_possorted_bam.md.bam \
    long_reads/GRCh38/NA24385_phased_possorted_bam.md.Hg38.bam
conda activate benchmarks

echo "$(tput setaf 1)$(tput setab 7)------- BAM files prepared (5.3/9.6) --------$(tput sgr 0)" 1>&3
