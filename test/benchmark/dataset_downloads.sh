# -------- -------- get data and unzip -------- -------- #
cd ../..
mkdir -p data && cd data

# mkdir -p short_reads && cd short_reads
# wget ...
# cd ..

echo "$(tput setaf 1)$(tput setab 7)------- short reads downloaded (3/7) --------$(tput sgr 0)" 1>&3

mkdir -p long_reads && cd long_reads
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam

wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam.bai
cd ..
echo "$(tput setaf 1)$(tput setab 7)------- long reads downloaded (4/7) --------$(tput sgr 0)" 1>&3

# -------- -------- get reference ../../data -------- -------- #
mkdir -p reference && cd reference

wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip --decompress --keep hs37d5.fa.gz
cd ..
echo "$(tput setaf 1)$(tput setab 7)------- reference downloaded (5/7) --------$(tput sgr 0)" 1>&3

samtools calmd -b data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam \
    > data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.md.bam data/reference/hs37d5.fa

echo "$(tput setaf 1)$(tput setab 7)------- missing MD tags added (6/7) --------$(tput sgr 0)" 1>&3

# -------- -------- get truth set ../../data -------- -------- #
mkdir -p truth_set && cd truth_set

wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/../../data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/../../data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/../../data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
cd ../..
echo "$(tput setaf 1)$(tput setab 7)------- truth set downloaded (7/7) --------$(tput sgr 0)" 1>&3
