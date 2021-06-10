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

for min_var_length in 10 30 50 100 500 1000 2000
do
    echo "Looping min_var_length number $min_var_length"
    /usr/bin/time -v ./build/iGenVar/bin/iGenVar \
        -j data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam \
        -o results/${current_date}/min_var_length_${min_var_length}_output.vcf -m 0 -m 1 \
        --min_var_length ${min_var_length}
done

echo "$(tput setaf 1)$(tput setab 7) run iGenVar with different min_var_length: 10 30 50 100 500 1000 2000 done (1/5) $(tput sgr 0)" 1>&3

for max_var_length in 500000 1000000 2000000
do
    echo "Looping max_var_length number $max_var_length"
    /usr/bin/time -v ./build/iGenVar/bin/iGenVar \
        -j data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam \
        -o results/${current_date}/max_var_length_${max_var_length}_output.vcf -m 0 -m 1 \
        --max_var_length ${max_var_length}
done

echo "$(tput setaf 1)$(tput setab 7) run iGenVar with different max_var_length: 500000 1000000 2000000 done (2/5) $(tput sgr 0)" 1>&3

for max_tol_inserted_length in 1 5 10
do
    echo "Looping max_tol_inserted_length number $max_tol_inserted_length"
    /usr/bin/time -v ./build/iGenVar/bin/iGenVar \
        -j data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam \
        -o results/${current_date}/max_tol_inserted_length_${max_tol_inserted_length}_output.vcf -m 0 -m 1 \
        --max_tol_inserted_length ${max_tol_inserted_length}
done

echo "$(tput setaf 1)$(tput setab 7) run iGenVar with different max_tol_inserted_length: 1 5 10 done (3/5) $(tput sgr 0)" 1>&3

for max_overlap in 5 10 20
do
    echo "Looping max_overlap number $max_overlap"
    /usr/bin/time -v ./build/iGenVar/bin/iGenVar \
        -j data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam \
        -o results/${current_date}/max_overlap_${max_overlap}_output.vcf -m 0 -m 1 \
        --max_overlap ${max_overlap}
done

echo "$(tput setaf 1)$(tput setab 7) run iGenVar with different max_overlap: 5 10 20 done (4/5) $(tput sgr 0)" 1>&3

for hierarchical_clustering_cutoff in 20 50 100 200 500 5000
do
    echo "Looping hierarchical_clustering_cutoff number $hierarchical_clustering_cutoff"
    /usr/bin/time -v ./build/iGenVar/bin/iGenVar \
        -j data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam \
        -o results/${current_date}/hc_cutoff_${hierarchical_clustering_cutoff}_output.vcf -m 0 -m 1 \
        --hierarchical_clustering_cutoff ${hierarchical_clustering_cutoff}
done

echo "$(tput setaf 1)$(tput setab 7) run iGenVar with different hierarchical_clustering_cutoff: 20 50 100 200 500 5000 done (5/5) $(tput sgr 0)" 1>&3

# grep VmPeak /proc/$PID/status

echo "$(tput setaf 1)$(tput setab 7)------- done: iGenVar with cigar string & split read method --------$(tput sgr 0)" 1>&3

end=$(date +%s)                 # get end-date
runtime=$(((end-start)/60))     # calculate runtime
echo "$runtime minutes"
echo "$(tput setaf 1)$(tput setab 7)This run took $runtime minutes$(tput sgr 0)" 1>&3
