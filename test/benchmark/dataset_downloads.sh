# -------- -------- get data and unzip -------- -------- #
mkdir -p short_reads && cd short_reads
# HG002 - NA24385:
# NIST Illumina 2x250bps
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam.bai

# NIST Stanford Illumina 6kb matepair
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG002.mate_pair.sorted.bam
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG002.mate_pair.sorted.bam.bai

cd ..

echo "$(tput setaf 1)$(tput setab 7)------- short reads downloaded (3.1/3.5) --------$(tput sgr 0)" 1>&3

mkdir -p long_reads && cd long_reads
# HG002 - NA24385:
## MtSinai PacBio (minimap2)
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh38.bam
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh38.bam.bai

## PacBio CCS 10kb (minimap2)
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam.bai

## 10X Genomics
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/10XGenomics/NA24385_phased_possorted_bam.bam
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/10XGenomics/NA24385_phased_possorted_bam.bam.bai

# HG005 - NA24631:
## MtSinai PacBio (minimap2)
# wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
#     https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/MtSinai_PacBio/PacBio_minimap2_bam/HG005_PacBio_GRCh38.bam

# wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
#     https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/MtSinai_PacBio/PacBio_minimap2_bam/HG005_PacBio_GRCh38.bam.bai

# ## PacBio SequelII CCS 11kb
# wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
#     https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/PacBio_SequelII_CCS_11kb/HG005.SequelII.pbmm2.hs37d5.whatshap.haplotag.10x.bam

# wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
#     https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/PacBio_SequelII_CCS_11kb/HG005.SequelII.pbmm2.hs37d5.whatshap.haplotag.10x.bam.bai

cd ..
echo "$(tput setaf 1)$(tput setab 7)------- long reads downloaded (3.2/3.5) --------$(tput sgr 0)" 1>&3

# -------- -------- get reference ../../data -------- -------- #
mkdir -p reference && cd reference

# Reference Genome (for MtSinai PacBio)
# GRCh38/hg38
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip --decompress --keep GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# 1000genomes Reference Genome Sequence (hs37d5), based on NCBI GRCh37 (for PacBio CCS 10kb)
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip --decompress --keep hs37d5.fa.gz

# ucsc (for 10X Genomics)
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
gzip --decompress --keep hg19.fa.gz
# reorder reference:
samtools faidx data/reference/hg19.fa $(cat Repos/iGenVar/test/benchmark/order_of_hg19.txt) > data/reference/hg19.reordered.fa

# GRCh37/hg19
#     ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
# GRCh38_reference_genome
#     https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd ..
echo "$(tput setaf 1)$(tput setab 7)------- reference downloaded (3.3/3.5) --------$(tput sgr 0)" 1>&3
# HG002 - NA24385:
## MtSinai PacBio (minimap2)
samtools calmd --threads 16 -Q -b data/long_reads/HG002_PacBio_GRCh38.bam \
    > data/long_reads/HG002_PacBio_GRCh38.md.bam \
    data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

## PacBio CCS 10kb
samtools calmd --threads 16 -Q -b data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam \
    > data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.md.bam \
    data/reference/GCA_000001405.28_GRCh38.p13_genomic.fna.gz

## 10X Genomics
samtools calmd -b data/long_reads/NA24385_phased_possorted_bam.bam \
    > data/long_reads/NA24385_phased_possorted_bam.md.bam \
    data/reference/hg38.fa

echo "$(tput setaf 1)$(tput setab 7)------- missing MD tags added (3.4/3.5) --------$(tput sgr 0)" 1>&3

# -------- -------- get truth set ../../data -------- -------- #
mkdir -p truth_set && cd truth_set

# HG002 - NA24385:
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed

# unzip vcf
bgzip -d -c data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz > data/truth_set/HG002_SVs_Tier1_v0.6.vcf
# add chr to chromosome names
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' data/truth_set/HG002_SVs_Tier1_v0.6.vcf > data/truth_set/HG002_SVs_Tier1_v0.6.renamed_chr.vcf
sed -e 's/^/chr/' data/truth_set/HG002_SVs_Tier1_v0.6.bed > data/truth_set/HG002_SVs_Tier1_v0.6.renamed_chr.bed
# zip file again:
bgzip -c data/truth_set/HG002_SVs_Tier1_v0.6.renamed_chr.vcf > data/truth_set/HG002_SVs_Tier1_v0.6.renamed_chr.vcf.gz
# create new tbi file
tabix -p vcf data/truth_set/HG002_SVs_Tier1_v0.6.renamed_chr.vcf.gz

# HG005 - NA24631: (there is no truth set so far)
# wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
#     https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/...     .vcf.gz
# wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
#     https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/...     .vcf.gz.tbi
# wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
#     https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/...     .bed

cd ..
echo "$(tput setaf 1)$(tput setab 7)------- truth set downloaded (3.5/3.5) --------$(tput sgr 0)" 1>&3
