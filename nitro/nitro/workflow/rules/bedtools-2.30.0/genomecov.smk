rule bedtools_genomecov_v2_30_0:
    input:
        bam = rules.samtools_merge_v1_15_1.output.merge,
        index = rules.samtools_index_v1_15_1.output.index,
    output:
        low_vcf = config["BaseDir"] / config["OutputDir"] / "output/{sample}/{sample}_no_cov.bed",
    params:
    threads: 1
    log:
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_no_cov.log"
    conda:
        "env.yaml"
    wrapper:
        "file:rules/bedtools-2.30.0/genomecov_wrapper.py"
