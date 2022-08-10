#!/usr/bin/env sh
echo "Logfile written to: initial_steps.log"

# -------- -------- Prepare logfile and runtime computation -------- -------- #
mkdir -p logs
mkdir -p logs/sniffles
mkdir -p logs/truvari
mkdir -p logs/caller_comparison_iGenVar_only

exec 3>&1 4>&2 # store original streams
trap 'exec 2>&4 1>&3' 0 1 2 3 # restore original streams when scripts ends
exec 1>logs/initial_steps.log 2>&1 # write stdout and stdcerr to logfile

set -ex

echo $(date)
start=$(date +%s) # get starting date

echo "$(tput setaf 1)$(tput setab 7)------- environment created (1/9) --------$(tput sgr 0)" 1>&3

# -------- -------- get & build iGenVar, SVIM, Vaquita, Vaquita-LR -------- -------- #
mkdir -p Repos && cd Repos/
git clone https://github.com/seqan/iGenVar.git
git clone https://github.com/eldariont/svim.git
git clone https://github.com/seqan/vaquita.git
git clone https://github.com/joshuak94/vaquita-LR.git # git@github.com:joshuak94/vaquita-LR.git

cd iGenVar
git submodule update --recursive --init

cd ../vaquita-LR
# Currently you have to update the .git/config git@ -> https://
git submodule update --recursive --init
git checkout develop

echo "$(tput setaf 1)$(tput setab 7)------- Repos downloaded (2/9) --------$(tput sgr 0)" 1>&3

# ToDo: Remove and use conda again when this is merged: https://github.com/bioconda/bioconda-recipes/pull/29099
# Install SVIM from github (requires Python 3.6.* or newer): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
cd ../svim
pip install .

echo "$(tput setaf 1)$(tput setab 7)------- SVIM installed (3/9) --------$(tput sgr 0)" 1>&3

cd ../../..
mkdir -p build && cd build && mkdir -p iGenVar && cd iGenVar
cmake ../../Repos/iGenVar/ -Weverything
make -j 16
make test && make doc

echo "$(tput setaf 1)$(tput setab 7)------- iGenVar built (4/9) --------$(tput sgr 0)" 1>&3

cd ..
mkdir -p Vaquita && cd Vaquita
cmake ../../Repos/vaquita/ -Weverything
make -j 16

echo "$(tput setaf 1)$(tput setab 7)------- Vaquita built (5/9) --------$(tput sgr 0)" 1>&3

cd ..
mkdir -p Vaquita-LR && cd Vaquita-LR
cmake ../../Repos/vaquita-LR/ -Weverything
make -j 16

echo "$(tput setaf 1)$(tput setab 7)------- Vaquita-LR built (6/9) --------$(tput sgr 0)" 1>&3

cd ../..
mkdir -p data && cd data

# short & long reads, reference and truth sets
./../Repos/iGenVar/test/benchmark/dataset_downloads.sh
# create simulated datasets
./../Repos/iGenVar/test/benchmark/simulation/mason_simulation.sh

echo "$(tput setaf 1)$(tput setab 7)------- Datasets downloaded (7/9) --------$(tput sgr 0)" 1>&3

# -------- -------- pre installation steps -------- -------- #
# We do our benchmarks with snakemake, and need several tools like
# python3, samtools, snakemake, tabix, pip and truvari
# we recommend to use anaconda: https://www.anaconda.com/products/individual
# and run
# wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
# ./Anaconda3-2021.05-Linux-x86_64.sh
conda env create -f Repos/iGenVar/test/benchmark/envs/environment.yml
conda activate benchmarks

echo "$(tput setaf 1)$(tput setab 7)------- anaconda environment prepared (8/9) --------$(tput sgr 0)" 1>&3

./../Repos/iGenVar/test/benchmark/prepare_BAM_with_crossmap.sh
./../Repos/iGenVar/test/benchmark/prepare_truth_set.sh
# ./../Repos/iGenVar/test/benchmark/prepare_truth_set_with_NCBI.sh

cd ..

echo "$(tput setaf 1)$(tput setab 7)------- data prepared (9/9) --------$(tput sgr 0)" 1>&3

mkdir -p tmp/picard
mkdir -p results

# run Snakefile with:
# snakemake --snakefile Repos/iGenVar/test/benchmark/caller_comparison_long_read/Snakefile --cores 16 --rerun-incomplete \
# --stats logs/{DATE}_snakemake_caller_comparison_benchmarks_stats.log \
# --runtime-profile logs/{DATE}_snake_caller_comparison_benchmarks_time.log --use-conda --conda-frontend mamba

# ----------------------------------------
echo "$(tput setaf 1)$(tput setab 7)------- Initial steps - done --------$(tput sgr 0)" 1>&3

end=$(date +%s)                                                         # get end-date
runtime=$(((end-start)/60))                                             # calculate runtime
echo "$runtime minutes"
echo "$(tput setaf 1)$(tput setab 7)This run took $runtime minutes$(tput sgr 0)" 1>&3
