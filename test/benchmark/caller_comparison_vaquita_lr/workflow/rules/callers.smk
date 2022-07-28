sample = config["parameters"]["sample"],
min_var_length = config["parameters"]["min_var_length"]

rule run_Vaquita_LR:
    output:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants_without_contigs.vcf"
    log:
        "logs/caller_comparison_vaquita_lr/{input_combination}_output.log"
    threads: 8
    run:
        # The combinations of S1L3, S2L1, S2L2 are not possible, as they have a different chromosome naming (with and without chr).
        if wildcards.input_combination == 'S1': # Illumina Paired End
            short_bam = config["short_read_bam"]["s1"],
            genome = config["reference_fa"]["Illumina_Paired_End"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --shortRead {short_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} &>> {log}
            """)
        elif wildcards.input_combination == 'S2': # Illumina Mate Pair
            short_bam = config["short_read_bam"]["s2"],
            genome = config["reference_fa"]["Illumina_Mate_Pair"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --shortRead {short_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} &>> {log}
            """)
        # [2022-7-27 14:13:41] [START] SViper: PREP PARALLEL PROCESSING
        # [2022-7-27 14:13:41] [END] SViper: PREP PARALLEL PROCESSING (0 seconds.)
        # [2022-7-27 14:13:41] [START] SViper: POLISHING VARIANTS
        # terminate called after throwing an instance of 'std::logic_error'
        # what():  Invalid region specified. Parameter regionStart was greater than regionEnd.
        # Command terminated by signal 6
        elif wildcards.input_combination == 'S1L1': # Illumina Paired End & MtSinai PacBio
            short_bam = config["short_read_bam"]["s1"],
            long_bam = config["long_read_bam"]["l1"],
            genome = config["reference_fa"]["Illumina_Paired_End"] # or config["reference_fa"]["MtSinai_PacBio"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --shortRead {short_bam} --longRead {long_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} --no-polishing &>> {log}
            """)
        # [2022-7-27 11:56:42] [START] SViper: PREP PARALLEL PROCESSING
        # [2022-7-27 11:56:42] [END] SViper: PREP PARALLEL PROCESSING (0 seconds.)
        # [2022-7-27 11:56:42] [START] SViper: POLISHING VARIANTS
        # Command terminated by signal 11
        elif wildcards.input_combination == 'S1L2': # Illumina Paired End & PacBio CCS
            short_bam = config["short_read_bam"]["s1"],
            long_bam = config["long_read_bam"]["l2"],
            genome = config["reference_fa"]["Illumina_Paired_End"] # or config["reference_fa"]["PacBio_CCS"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --shortRead {short_bam} --longRead {long_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} --no-polishing &>> {log}
            """)
        elif wildcards.input_combination == 'S2L3': # Illumina Mate Pair & 10X Genomics
            short_bam = config["short_read_bam"]["s2"],
            long_bam = config["long_read_bam"]["l3"],
            genome = config["reference_fa"]["Illumina_Mate_Pair"] # or config["reference_fa"]["10X_Genomics"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --shortRead {short_bam} --longRead {long_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} &>> {log}
            """)
        elif wildcards.input_combination == 'L1': # MtSinai PacBio
            long_bam = config["long_read_bam"]["l1"],
            genome = config["reference_fa"]["MtSinai_PacBio"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --longRead {long_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} &>> {log}
            """)
        elif wildcards.input_combination == 'L2': # PacBio CCS
            long_bam = config["long_read_bam"]["l2"],
            genome = config["reference_fa"]["PacBio_CCS"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --longRead {long_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} &>> {log}
            """)
        else: # wildcards.input_combination == 'L3': # 10X Genomics
            long_bam = config["long_read_bam"]["l3"],
            genome = config["reference_fa"]["10X_Genomics"]
            shell("""
            /usr/bin/time -v ./build/Vaquita-LR/bin/vaquita call --longRead {long_bam} \
                --referenceGenome {genome} --cutoff 1 --minSVSize {min_var_length} --threadCount {threads} \
                --output-file {output.vcf} &>> {log}
            """)
        # Defaults:
        # GENERAL
        #     -v, --minVote (signed 32 bit integer)
        #         Minimum number of evidence types(=vote) that support SVs for rescue. -1: Supported by all evidence types.
        #         Default: -1. Value must be in range [-1,2.14748e+09].
        #     -q, --minMapQual (signed 32 bit integer)
        #         Mapping quaility cutoff. Default: 20. Value must be in range [0,2.14748e+09].
        #     -a, --adjTol (signed 32 bit integer)
        #         Positional adjacency in nucleotide resolution. Default: 50. Value must be in range [0,2.14748e+09].
        #     --no-polishing
        #         Do not polish the variants based on SViper. [Only relevant when short and long reads are both provided.]
        #     -w, --write-breakpoint
        #         Write breakpoint informtation in a tab-separated format (breakpoints.tsv).
        # SPLIT-READ EVIDENCE
        #     -s, --maxSplit (signed 32 bit integer)
        #         Maximum number of segments in a single read. Default: 2. Value must be in range [0,2.14748e+09].
        #     --maxSplitLong (signed 32 bit integer)
        #         Maximum number of segments in a single long read. Only relevant for long read input (-z) Default: 7. Value
        #         must be in range [0,2.14748e+09].
        #     -x, --maxOverlap (signed 32 bit integer)
        #         Maximum allowed overlaps between segments. Default: 20. Value must be in range [0,2.14748e+09].
        #     --maxOverlapLong (signed 32 bit integer)
        #         Maximum allowed overlaps between long read segments. Only relevant for long read input (-z) Default: 1000.
        #         Value must be in range [0,2.14748e+09].
        #     -b, --minSplitReadSupport (double)
        #         SVs supported by >= b get a vote. Default: 1. Value must be in range [0,1.79769e+308].
        # READ-PAIR EVIDENCE
        #     -p, --pairedEndSearchSize (signed 32 bit integer)
        #         Size of breakpoint candidate regions. Default: 500. Value must be in range [0,2.14748e+09].
        #     -i, --abInsParam (double)
        #         Discordant insertion size: median +/- (MAD * i) Default: 9. Value must be in range [0,1.79769e+308].
        #     -d, --depthOutlier (double)
        #         Depth outlier: {Q3 + (IQR * d)} Default: 1. Value must be in range [0,1.79769e+308].
        #     -e, --minPairSupport (double)
        #         SVs supported by >= e get a vote. Default: 1. Value must be in range [0,1.79769e+308].
        # SOFT-CLIPPED EVIDENCE
        #     -l, --minClippedSeqSize (signed 32 bit integer)
        #         Minimum size of clipped sequence to be considered. Default: 20. Value must be in range [0,2.14748e+09].
        #     -t, --clippedSeqErrorRate (double)
        #         Maximum edit distance: floor{length of clipped sequence * (1 - t)}. Default: 0.1. Value must be in range
        #         [0,1.79769e+308].
        # READ-DEPTH EVIDENCE
        #     -n, --samplingNum (signed 32 bit integer)
        #         Number of random sample to estimate the background distribution(Q3, IQR, ..) of read-depth evidence.
        #         Default: 100000. Value must be in range [0,2.14748e+09].
        #     -f, --readDepthWindowSize (signed 32 bit integer)
        #         Window size to caclulate average read-depth around breakpoints. Default: 20. Value must be in range
        #         [0,2.14748e+09].
        #     --use-re-for-bs
        #         Use RE for balanced SVs(eg. inverison).
        #     -j, --reThreshold (double)
        #         SVs satisfy read-depth evidence >= {Q3 + (IQR * h)} get a vote. Default: 1. Value must be in range
        #         [0,1.79769e+308].

rule fix_vaquita_lr:
    input:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants_without_contigs.vcf"
    output:
        vcf_1 = "results/caller_comparison_vaquita_lr/{input_combination}/variants_without_fileformat.vcf",
        vcf_2 = "results/caller_comparison_vaquita_lr/{input_combination}/variants.vcf"
    params:
        igenvar_vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.vcf"
    run:
        shell("""
        fileformat=$(head -n 1 {input.vcf}) && \
        less {params.igenvar_vcf} | grep contig | cat - {input.vcf} > {output.vcf_1} && \
        echo $fileformat | cat - {output.vcf_1} > {output.vcf_2}
        """)
        # without the 'chr' prefix: S1, L1, L2, S1L1, S1L2
        if (wildcards.input_combination == 'S1') | (wildcards.input_combination == 'L1') | (wildcards.input_combination == 'L2') | (wildcards.input_combination == 'S1L1') | (wildcards.input_combination == 'S1L2'):
            shell("sed -i -e 's/ID=chr/ID=/g' {output.vcf_2}")
