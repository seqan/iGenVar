min_qual_iGenVar = list(range(config["quality_ranges"]["iGenVar"]["from"],
                              config["quality_ranges"]["iGenVar"]["to"],
                              config["quality_ranges"]["iGenVar"]["step"]))

rule filter_vcf:
    input:
        "results/caller_comparison_short_read/{caller,iGenVar_S|iGenVar_SL}/variants.vcf"
    output:
        "results/caller_comparison_short_read/{caller,iGenVar_S|iGenVar_SL}/variants.min_qual_{min_qual}.vcf"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}' {input} > {output}"

rule bgzip:
    input:
        "{name}.vcf"
    output:
        "{name}.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule truvari:
    input:
        vcf = "results/caller_comparison_short_read/{caller}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_short_read/{caller}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    params:
        output_dir = "results/caller_comparison_short_read/eval/{caller}/min_qual_{min_qual}",
        truth_set_gz = config["truth_set"]["gz"],
        truth_set_bed = config["truth_set"]["bed"]
    output:
        summary = "results/caller_comparison_short_read/eval/{caller}/min_qual_{min_qual}/summary.txt"
    shell:
        """
        rm -rf {params.output_dir} && truvari bench -b {params.truth_set_gz} -c {input.vcf} -o {params.output_dir} \
            --passonly --includebed {params.truth_set_bed} -p 0 &>> logs/truvari_output.log
        """

rule reformat_truvari_results:
    input:
        "results/caller_comparison_short_read/eval/{caller}/min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_short_read/eval/{caller}/min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule cat_truvari_results_all:
    input:
        iGenVar_S   = expand("results/caller_comparison_short_read/eval/iGenVar_S/min_qual_{min_qual}/pr_rec.txt",
                             min_qual = min_qual_iGenVar),
        iGenVar_SL  = expand("results/caller_comparison_short_read/eval/iGenVar_SL/min_qual_{min_qual}/pr_rec.txt",
                             min_qual = min_qual_iGenVar)
    output:
        iGenVar_S   = temp("results/caller_comparison_short_read/eval/igenvar_s.all_results.txt"),
        iGenVar_SL  = temp("results/caller_comparison_short_read/eval/igenvar_sl.all_results.txt")
        all = "results/caller_comparison_short_read/eval/all_results.txt"
    threads: 1
    run:
        shell("cat {input.iGenVar_S} > {output.iGenVar_S}")
        shell("cat {input.iGenVar_SL} > {output.iGenVar_SL}")
        shell("cat {output.iGenVar_S} {output.iGenVar_SL} > {output.all}")
