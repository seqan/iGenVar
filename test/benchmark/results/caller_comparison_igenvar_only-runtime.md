# 2022-03-30 & 2022-04-20

## L1

Done with clustering. Found 8770966 junction clusters.
No refinement was selected.
Detected 4816877 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_long_reads data/long_reads/GRCh37/HG002_PacBio_GRCh37.bam --output results/caller_comparison_iGenVar_only/L1/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 7861.94
        System time (seconds): 2660.63
        Percent of CPU this job got: 213%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:22:15
        Maximum resident set size (kbytes): 19226000
        Major (requiring I/O) page faults: 558
        Minor (reclaiming a frame) page faults: 29390380
        Voluntary context switches: 403842
        Involuntary context switches: 1630405
        Swaps: 0
        File system inputs: 196703192
        File system outputs: 6534592

## L2

Done with clustering. Found 67693 junction clusters.
No refinement was selected.
Detected 67693 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_long_reads data/long_reads/GRCh37/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam --output results/caller_comparison_iGenVar_only/L2/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 2894.02
        System time (seconds): 782.43
        Percent of CPU this job got: 156% - 175%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 38:29.40 - 39:14.92
        Maximum resident set size (kbytes): 760840 - 761000
        Major (requiring I/O) page faults: 2
        Minor (reclaiming a frame) page faults: 483754
        Voluntary context switches: 230340
        Involuntary context switches: 282670339
        Swaps: 0
        File system inputs: 121766248
        File system outputs: 1293312

## L3

Done with clustering. Found 6314105 junction clusters.
No refinement was selected.
Detected 368354 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_long_reads data/long_reads/GRCh37/NA24385_phased_possorted_bam.bam --output results/caller_comparison_iGenVar_only/L3/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 15293.00
        System time (seconds): 8260.35
        Percent of CPU this job got: 204% - 249%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:05:09 - 3:12:03
        Maximum resident set size (kbytes): 3489880 - 3489888
        Major (requiring I/O) page faults: 7113
        Minor (reclaiming a frame) page faults: 3178988
        Voluntary context switches: 191910
        Involuntary context switches: 1901718237
        Swaps: 0
        File system inputs: 105495784
        File system outputs: 95754136

## S1

Done with clustering. Found 29029 junction clusters.
No refinement was selected.
Detected 29029 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.hs37d5.2x250.sorted.bam --output results/caller_comparison_iGenVar_only/S1/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 11302.52
        System time (seconds): 5329.20
        Percent of CPU this job got: 261%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:46:10
        Maximum resident set size (kbytes): 855960
        Major (requiring I/O) page faults: 505
        Minor (reclaiming a frame) page faults: 322586
        Voluntary context switches: 167182
        Involuntary context switches: 1692902
        Swaps: 0
        File system inputs: 89807880
        File system outputs: 174112536

## S2

Done with clustering. Found 1158118 junction clusters.
No refinement was selected.
Detected 7952 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.mate_pair.sorted.bam --output results/caller_comparison_iGenVar_only/S2/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 5925.02
        System time (seconds): 3763.05
        Percent of CPU this job got: 264% - 297%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:01:01 - 2:00:45
        Maximum resident set size (kbytes): 2814744 - 3423272
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 24
        Minor (reclaiming a frame) page faults: 1450353
        Voluntary context switches: 56288
        Involuntary context switches: 41373403
        Swaps: 0
        File system inputs: 33301568
        File system outputs: 83639696

## S1L1

Done with clustering. Found 8779621 junction clusters.
No refinement was selected.
Detected 4825813 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.hs37d5.2x250.sorted.bam --input_long_reads data/long_reads/GRCh37/HG002_PacBio_GRCh37.bam --output results/caller_comparison_iGenVar_only/S1L1/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 34941.41
        System time (seconds): 15371.73
        Percent of CPU this job got: 336%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 4:09:12
        Maximum resident set size (kbytes): 19296920
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 27259
        Minor (reclaiming a frame) page faults: 45705194
        Voluntary context switches: 744697
        Involuntary context switches: 59971300
        Swaps: 0
        File system inputs: 331492688
        File system outputs: 181245928

## S1L2

Done with clustering. Found 87137 junction clusters.
No refinement was selected.
Detected 87137 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.hs37d5.2x250.sorted.bam --input_long_reads data/long_reads/GRCh37/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam --output results/caller_comparison_iGenVar_only/S1L2/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 14321.41
        System time (seconds): 6053.73
        Percent of CPU this job got: 265%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:08:02
        Maximum resident set size (kbytes): 973364
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 1943
        Minor (reclaiming a frame) page faults: 513648
        Voluntary context switches: 231279
        Involuntary context switches: 1915351
        Swaps: 0
        File system inputs: 109554536
        File system outputs: 175416232

## S1L3

Done with clustering. Found 6342150 junction clusters.
No refinement was selected.
Detected 396401 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.hs37d5.2x250.sorted.bam --input_long_reads data/long_reads/GRCh37/NA24385_phased_possorted_bam.bam --output results/caller_comparison_iGenVar_only/S1L3/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 21123.33
        System time (seconds): 9437.12
        Percent of CPU this job got: 271%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 3:07:19
        Maximum resident set size (kbytes): 3567316
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 122
        Minor (reclaiming a frame) page faults: 2009117
        Voluntary context switches: 216085
        Involuntary context switches: 2784023
        Swaps: 0
        File system inputs: 130026216
        File system outputs: 269948488

## S2L1

Done with clustering. Found 9925351 junction clusters.
No refinement was selected.
Detected 4822602 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.mate_pair.sorted.bam --input_long_reads data/long_reads/GRCh37/HG002_PacBio_GRCh37.bam --output results/caller_comparison_iGenVar_only/S2L1/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 13539.96
        System time (seconds): 6196.35
        Percent of CPU this job got: 249%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:11:42
        Maximum resident set size (kbytes): 19670912
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 4169
        Minor (reclaiming a frame) page faults: 20250986
        Voluntary context switches: 425718
        Involuntary context switches: 43271983
        Swaps: 0
        File system inputs: 208965424
        File system outputs: 90174440

## S2L2

Done with clustering. Found 1224088 junction clusters.
No refinement was selected.
Detected 73938 SVs.
        Command being timed: `"./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.mate_pair.sorted.bam --input_long_reads data/long_reads/GRCh37/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam --output results/caller_comparison_iGenVar_only/S2L2/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 15647.20
        System time (seconds): 9037.36
        Percent of CPU this job got: 248% - 253%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:28:24 - 2:42:28
        Maximum resident set size (kbytes): 2828092 - 3417652
        Major (requiring I/O) page faults: 42
        Minor (reclaiming a frame) page faults: 10475699
        Voluntary context switches: 389753
        Involuntary context switches: 1398807520
        Swaps: 0
        File system inputs: 192060304
        File system outputs: 86744600

## S2L3

Done with clustering. Found 7457031 junction clusters.
No refinement was selected.
Detected 373861 SVs.
        Command being timed: `./build/iGenVar/bin/iGenVar --input_short_reads data/short_reads/GRCh37/HG002.mate_pair.sorted.bam --input_long_reads data/long_reads/GRCh37/NA24385_phased_possorted_bam.bam --output results/caller_comparison_iGenVar_only/S2L3/variants.vcf --vcf_sample_name HG002 --threads 2 --verbose --min_var_length 40 --max_var_length 1000000 --min_qual 1`
        User time (seconds): 25386.84
        System time (seconds): 14474.91
        Percent of CPU this job got: 253% - 311%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:49:55 - 4:21:37
        Maximum resident set size (kbytes): 4173176 - 4173960
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 37214
        Minor (reclaiming a frame) page faults: 5393501
        Voluntary context switches: 234848
        Involuntary context switches: 2703304840
        Swaps: 0
        File system inputs: 104955264
        File system outputs: 180960336
