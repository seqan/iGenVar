rule filter_SVIM_vcf:
    input:
        vcf = "results/caller_comparison_long_read/{dataset}/{caller}/variants.vcf"
    output:
        vcfs = "results/caller_comparison_long_read/{dataset}/{caller}/variants.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}' {input.vcf} > {output.vcfs}"

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

# do not use from bioconda: truvari 0.1.2018.08.10 (python 2.7)
# use from pip: Truvari v3.0.0 (python 3.6)
# The input .gz file needs to be compressed by bgzip.
rule truvari:
    input:
        vcf = "results/caller_comparison_long_read/{dataset}/{caller}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_long_read/{dataset}/{caller}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_long_read/{dataset}/eval/{caller}/min_qual_{min_qual}/summary.txt"
    log:
        "logs/caller_comparison_long_read/truvari/truvari_output.{dataset}.{caller}.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_long_read/{dataset}/eval/{caller}/min_qual_{min_qual}"
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            truth_set_gz = config["truth_set"]["gz"],
            truth_set_bed = config["truth_set"]["bed"]
        elif wildcards.dataset == 'PacBio_CCS':
            truth_set_gz = config["truth_set"]["gz"],
            truth_set_bed = config["truth_set"]["bed"]
        else: # wildcards.dataset == '10X_Genomics'
            truth_set_gz = config["truth_set_renamed_chr"]["gz"],
            truth_set_bed = config["truth_set_renamed_chr"]["bed"]
        shell("""
            rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
                --passonly --includebed {truth_set_bed} &>> {log}
        """)
        # -f data/reference/hs37d5.fa

rule reformat_truvari_results:
    input:
        "results/caller_comparison_long_read/{dataset}/eval/{caller}/min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_long_read/{dataset}/eval/{caller}/min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.min_qual}\", $1, $2 }}' > {output}
        """

rule cat_truvari_results_all:
    input:
        igenvar_L        = expand("results/caller_comparison_long_read/{{dataset}}/eval/iGenVar_L/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                                      config["quality_ranges"]["iGenVar"]["to"],
                                                      config["quality_ranges"]["iGenVar"]["step"]))),
        igenvar_SL       = expand("results/caller_comparison_long_read/{{dataset}}/eval/iGenVar_SL/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                                      config["quality_ranges"]["iGenVar"]["to"],
                                                      config["quality_ranges"]["iGenVar"]["step"]))),
        svim             = expand("results/caller_comparison_long_read/{{dataset}}/eval/SVIM/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"],
                                                      config["quality_ranges"]["step"]))),
        sniffles         = expand("results/caller_comparison_long_read/{{dataset}}/eval/Sniffles/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["sniffles"]["from"],
                                                      config["quality_ranges"]["sniffles"]["to"],
                                                      config["quality_ranges"]["sniffles"]["step"]))),
        pbsv             = expand("results/caller_comparison_long_read/{{dataset}}/eval/pbsv/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"],
                                                      config["quality_ranges"]["step"]))),
        pbsv_without_DUP = expand("results/caller_comparison_long_read/{{dataset}}/eval/pbsv_without_DUP/min_qual_{min_qual}/pr_rec.txt",
                                  min_qual=list(range(config["quality_ranges"]["from"],
                                                      config["quality_ranges"]["to"],
                                                      config["quality_ranges"]["step"])))
    output:
        igenvar_L        = temp("results/caller_comparison_long_read/{{dataset}}/eval/igenvar_L.all_results.txt"),
        igenvar_SL       = temp("results/caller_comparison_long_read/{{dataset}}/eval/igenvar_SL.all_results.txt"),
        svim             = temp("results/caller_comparison_long_read/{{dataset}}/eval/svim.all_results.txt"),
        sniffles         = temp("results/caller_comparison_long_read/{{dataset}}/eval/sniffles.all_results.txt"),
        pbsv             = temp("results/caller_comparison_long_read/{{dataset}}/eval/pbsv.all_results.txt"),
        pbsv_without_DUP = temp("results/caller_comparison_long_read/{{dataset}}/eval/pbsv_without_DUP.all_results.txt"),
        all = "results/caller_comparison_long_read/{dataset}/eval/all_results.txt"
    threads: 1
    run:
        shell("cat {input.igenvar_L} > {output.igenvar_L}")
        shell("cat {input.igenvar_SL} > {output.igenvar_SL}")
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {input.pbsv_without_DUP} > {output.pbsv_without_DUP}")
        shell("cat {output.igenvar_L} {output.igenvar_SL} {output.svim} \
                   {output.sniffles} {output.pbsv} {output.pbsv_without_DUP} > {output.all}")
