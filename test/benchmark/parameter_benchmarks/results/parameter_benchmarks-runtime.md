# Parameter benchmarks

# 2021-09-17

"total_runtime": 771.5161242485046,
"rules": {
    "run_igenvar": {
        "mean-runtime": 690.0873486995697,
        "min-runtime": 671.3452603816986,
        "max-runtime": 708.8294370174408
    },

name                                        ncall       tsub      ttot      tavg
... /parameter_benchmarks/Snakefile:24 __rule_run_igenvar                2           0.000261  0.019935  0.009967

# 2021-09-15

"total_runtime": 3300.3677134513855,
"rules": {
    "run_igenvar": {
        "mean-runtime": 2468.2231792041234,
        "min-runtime": 1933.8672742843628,
        "max-runtime": 2733.89550113678
    },

name                                        ncall       tsub      ttot      tavg
... /parameter_benchmarks/Snakefile:20 __rule_run_igenvar                21          0.001778  0.166659  0.007936

## 2021-06-08
This run took 270 minutes

## 2021-06-10
This run took 271 minutes

    Command being timed: "./build/iGenVar/bin/iGenVar -j data/long_reads/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio_sorted.bam -o results/2021-06-10/hc_cutoff_5000_output.vcf -m 0 -m 1 --hierarchical_clustering_cutoff 5000"
    User time (seconds): 33649.52
    System time (seconds): 5535.41
    Percent of CPU this job got: 2839%
    Elapsed (wall clock) time (h:mm:ss or m:ss): 23:00.04
    Average shared text size (kbytes): 0
    Average unshared data size (kbytes): 0
    Average stack size (kbytes): 0
    Average total size (kbytes): 0
    Maximum resident set size (kbytes): 2112596
    Average resident set size (kbytes): 0
    Major (requiring I/O) page faults: 0
    Minor (reclaiming a frame) page faults: 278024
    Voluntary context switches: 1137
    Involuntary context switches: 291007535
    Swaps: 0
    File system inputs: 0
    File system outputs: 219936
    Socket messages sent: 0
    Socket messages received: 0
    Signals delivered: 0
    Page size (bytes): 4096
    Exit status: 0
