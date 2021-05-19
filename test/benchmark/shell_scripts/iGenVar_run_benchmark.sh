#!/usr/bin/env sh
current_date=$(date +"%Y-%m-%d")
echo "Logifle written to: ${current_date}_run_benchmark.log"

# -------- -------- Prepare logfile and runtime computation -------- -------- #

exec 3>&1 4>&2 # store original streams
trap 'exec 2>&4 1>&3' 0 1 2 3 # restore original streams when scripts ends
exec 1>logs/${current_date}_run_benchmark.log 2>&1 # write stdout and stdcerr to logfile

set -ex

echo $(date)
start=$(date +%s) # get starting date

# -------- -------- run iGenVar with default values -------- -------- #

echo "$(tput setaf 1)$(tput setab 7)------- run iGenVar with Methods: 1, 2, 3 --------$(tput sgr 0)" 1>&3

mkdir -p results/${current_date}

for min_qual in 1 5 10
do
    echo "Looping min_qual number $min_qual"

    for hierarchical_clustering_cutoff in 10 50 100 150
    do
        echo "Looping hierarchical_clustering_cutoff number $hierarchical_clustering_cutoff"
        /usr/bin/time -v ./build/iGenVar/bin/iGenVar \
            -j data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam \
            -o results/${current_date}/${min_qual}_${hierarchical_clustering_cutoff}_output.vcf -m 0 -m 1 \
            --min_qual ${min_qual} --hierarchical_clustering_cutoff ${hierarchical_clustering_cutoff}
    done
done

# grep VmPeak /proc/$PID/status

echo "$(tput setaf 1)$(tput setab 7)------- done: iGenVar with cigar string & split read method --------$(tput sgr 0)" 1>&3

end=$(date +%s)                 # get end-date
runtime=$(((end-start)/60))     # calculate runtime
echo "$runtime minutes"
echo "$(tput setaf 1)$(tput setab 7)This run took $runtime minutes$(tput sgr 0)" 1>&3
