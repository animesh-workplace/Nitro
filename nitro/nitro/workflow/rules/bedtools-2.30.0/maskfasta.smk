rule bedtools_maskfasta_v2_30_0:
    input:
        no_vcf = rules.bedtools_genomecov_v2_30_0.output.low_vcf,
    params:
        reference = "resources/reference/NC_045512.2.fna",
    output:
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/{sample}.reference.masked.fasta",
    threads: 1
    conda: "env.yaml"
    log:
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_consensus_masked.log"
    wrapper:
        "file:rules/bedtools-2.30.0/maskfasta_wrapper.py"
