#!/usr/bin/env sh
echo "Logifle written to: initial_steps.log"

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

cd iGenVar
git submodule update --recursive --init

echo "$(tput setaf 1)$(tput setab 7)------- iGenVar downloaded (1/4) --------$(tput sgr 0)" 1>&3

cd ../..
mkdir -p build && cd build && mkdir -p iGenVar && cd iGenVar
cmake ../../Repos/iGenVar/ -Weverything

make -j 16

make test && make doc

echo "$(tput setaf 1)$(tput setab 7)------- iGenVar built (2/4) --------$(tput sgr 0)" 1>&3

# -------- -------- get data and unzip -------- -------- #
cd ../..
mkdir -p data && cd data

# mkdir -p short_reads && cd short_reads
# wget ...
# cd ..

echo "$(tput setaf 1)$(tput setab 7)------- short reads downloaded (3/4) --------$(tput sgr 0)" 1>&3

mkdir -p long_reads && cd long_reads
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam

echo "$(tput setaf 1)$(tput setab 7)------- long reads downloaded (4/4) --------$(tput sgr 0)" 1>&3

cd ../..
mkdir -p results

# ----------------------------------------
echo "$(tput setaf 1)$(tput setab 7)------- Initial steps - done --------$(tput sgr 0)" 1>&3

end=$(date +%s)                                                         # get end-date
runtime=$(((end-start)/60))                                             # calculate runtime
echo "$runtime minutes"
echo "$(tput setaf 1)$(tput setab 7)This run took $runtime minutes$(tput sgr 0)" 1>&3
