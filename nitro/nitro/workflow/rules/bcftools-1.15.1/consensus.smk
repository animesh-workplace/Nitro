rule bcftools_consensus_v1_15_1:
    input:
        vcf = rules.gatk_haplotype_v4_2_6_1.output,
        mask = rules.bedtools_genomecov_v2_30_0.output.low_vcf
    params:
        reference = "resources/reference/NC_045512.2.fna",
    output:
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/{sample}.consensus.fasta",
    threads: 1
    conda: "env.yaml"
    log:
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_consensus.log"
    wrapper:
        "file:rules/bcftools-1.15.1/consensus_wrapper.py"
