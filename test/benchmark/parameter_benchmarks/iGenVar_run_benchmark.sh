#!/usr/bin/env sh
current_date=$(date +"%Y-%m-%d")
echo "Logfile written to: ${current_date}_run_benchmark.log"

# -------- -------- Prepare logfile and runtime computation -------- -------- #

exec 3>&1 4>&2 # store original streams
trap 'exec 2>&4 1>&3' 0 1 2 3 # restore original streams when scripts ends
exec 1>logs/${current_date}_run_benchmark.log 2>&1 # write stdout and stdcerr to logfile

set -ex

echo $(date)
start=$(date +%s) # get starting date

# -------- -------- pre installation steps -------- -------- #
# We do our benchmarks with snakemake, and need
# python3, samtools, snakemake, tabix, pip and truvari
# we recommend to use miniconda: https://docs.conda.io/en/latest/miniconda.html
# and run
# -------- -------- --------  -------- -------- -------- #
# conda env create -f environment.yml

# source ~/.bashrc                  # if the .bashrc is not executed automatically
conda activate benchmarks

# -------- -------- run iGenVar with default values -------- -------- #

echo "$(tput setaf 1)$(tput setab 7)------- run Snakefile --------$(tput sgr 0)" 1>&3

cores=1

if [ $1 == "--cores" ]
then
   echo "Use $2 cores."
   a = $2
else
   echo "Use one core."
fi

/usr/bin/time -v snakemake --snakefile Repos/iGenVar/test/benchmark/parameter_benchmarks/Snakefile --cores $cores \ # --quiet --forcerun run_igenvar \
--rerun-incomplete --stats logs/${current_date}_snakemake_stats.txt --timestamp --use-conda

echo "$(tput setaf 1)$(tput setab 7)------- done: iGenVar with cigar string & split read method --------$(tput sgr 0)" 1>&3

end=$(date +%s)                 # get end-date
runtime=$(((end-start)/60))     # calculate runtime
echo "$runtime minutes"
echo "$(tput setaf 1)$(tput setab 7)This run took $runtime minutes$(tput sgr 0)" 1>&3
