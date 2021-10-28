rule filter_vcf:
    input:
        "results/caller_comparison/{caller,iGenVar|SVIM}/variants.vcf"
    output:
        "results/caller_comparison/{caller,iGenVar|SVIM}/variants.min_qual_{min_qual}.vcf"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}' {input} > {output}"

rule bgzip:
    input:
        "{name}.vcf"
    output:
        "{name}.vcf.gz"
    shell:
        "pbgzip -c {input} > {output}"

rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule truvari:
    input:
        vcf = "results/caller_comparison/{caller}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison/{caller}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    params:
        output_dir = "results/caller_comparison/eval/{caller}/min_qual_{min_qual}"
    output:
        summary = "results/caller_comparison/eval/{caller}/min_qual_{min_qual}/summary.txt"
    shell:
        """
        rm -rf {params.output_dir} && truvari bench -b data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz \
        -c {input.vcf} -o {params.output_dir} --passonly --includebed data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0 \
        &>> logs/truvari_output.log
        """

rule reformat_truvari_results:
    input:
        "results/caller_comparison/eval/{caller}/min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison/eval/{caller}/min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
        | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.min_qual}\", $1, $2 }}' > {output}
        """

rule cat_truvari_results_all:
    input:
        igenvar          = expand("results/caller_comparison/eval/iGenVar/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"]+1,
                                                      config["quality_ranges"]["step"]))),
        svim             = expand("results/caller_comparison/eval/SVIM/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"]+1,
                                                      config["quality_ranges"]["step"]))),
        sniffles         = expand("results/caller_comparison/eval/Sniffles/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"]+1,
                                                      config["quality_ranges"]["step"]))),
        pbsv             = expand("results/caller_comparison/eval/pbsv/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"]+1,
                                                      config["quality_ranges"]["step"]))),
        pbsv_without_DUP = expand("results/caller_comparison/eval/pbsv_without_DUP/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"]+1,
                                                      config["quality_ranges"]["step"])))
    output:
        igenvar          = temp("results/caller_comparison/eval/igenvar.all_results.txt"),
        svim             = temp("results/caller_comparison/eval/svim.all_results.txt"),
        sniffles         = temp("results/caller_comparison/eval/sniffles.all_results.txt"),
        pbsv             = temp("results/caller_comparison/eval/pbsv.all_results.txt"),
        pbsv_without_DUP = temp("results/caller_comparison/eval/pbsv_without_DUP.all_results.txt"),
        all = "results/caller_comparison/eval/all_results.txt"
    threads: 1
    run:
        shell("cat {input.igenvar} > {output.igenvar}")
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {input.pbsv_without_DUP} > {output.pbsv_without_DUP}")
        shell("cat {output.igenvar} {output.svim} {output.sniffles} {output.pbsv} {output.pbsv_without_DUP} > {output.all}")
