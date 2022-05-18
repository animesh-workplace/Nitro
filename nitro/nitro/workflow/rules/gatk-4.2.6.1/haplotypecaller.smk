rule gatk_haplotype_v4_2_6_1:
    input:
        bam = rules.samtools_merge_v1_15_1.output.merge,
        index = rules.samtools_index_v1_15_1.output.index,
        # reference = rules.bedtools_maskfasta_v2_30_0.output,
    output:
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/{sample}.vcf.gz",
    params:
        reference = "resources/reference/NC_045512.2.fna"
    threads: 1
    log: 
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_haplotypecaller.log"
    conda:
        "env.yaml"
    wrapper:
        "file:rules/gatk-4.2.6.1/haplotypecaller_wrapper.py"
