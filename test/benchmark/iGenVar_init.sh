#!/usr/bin/env sh
echo "Logfile written to: initial_steps.log"

# -------- -------- Prepare logfile and runtime computation -------- -------- #
mkdir -p logs

exec 3>&1 4>&2 # store original streams
trap 'exec 2>&4 1>&3' 0 1 2 3 # restore original streams when scripts ends
exec 1>logs/initial_steps.log 2>&1 # write stdout and stdcerr to logfile

set -ex

echo $(date)
start=$(date +%s) # get starting date

# -------- -------- get & build iGenVar -------- -------- #
mkdir -p Repos && cd Repos/
git clone https://github.com/seqan/iGenVar.git
git clone https://github.com/eldariont/svim.git

cd iGenVar
git submodule update --recursive --init

# ToDo: Remove and use conda again when this is merged: https://github.com/bioconda/bioconda-recipes/pull/29099
# Install SVIM from github (requires Python 3.6.* or newer): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
cd ../svim
pip install .

echo "$(tput setaf 1)$(tput setab 7)------- iGenVar downloaded (1/7) --------$(tput sgr 0)" 1>&3

cd ../..
mkdir -p build && cd build && mkdir -p iGenVar && cd iGenVar
cmake ../../Repos/iGenVar/ -Weverything

make -j 16

make test && make doc

echo "$(tput setaf 1)$(tput setab 7)------- iGenVar built (2/7) --------$(tput sgr 0)" 1>&3

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

mkdir -p results

# -------- -------- pre installation steps -------- -------- #
# We do our benchmarks with snakemake, and need several tools like
# python3, samtools, snakemake, tabix, pip and truvari
# we recommend to use anaconda: https://www.anaconda.com/products/individual
# and run
# wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
# ./Anaconda3-2021.05-Linux-x86_64.sh
conda env create -f Repos/iGenVar/test/benchmark/envs/environment.yml

# run Snakefile with:
# snakemake --snakefile Repos/iGenVar/test/benchmark/caller_comparison/Snakefile --cores 16 --rerun-incomplete \
# --stats logs/{DATE}_snakemake_caller_comparison_benchmarks_stats.log \
# --runtime-profile logs/{DATE}_snake_caller_comparison_benchmarks_time.log --use-conda --conda-frontend mamba

# ----------------------------------------
echo "$(tput setaf 1)$(tput setab 7)------- Initial steps - done --------$(tput sgr 0)" 1>&3

end=$(date +%s)                                                         # get end-date
runtime=$(((end-start)/60))                                             # calculate runtime
echo "$runtime minutes"
echo "$(tput setaf 1)$(tput setab 7)This run took $runtime minutes$(tput sgr 0)" 1>&3
