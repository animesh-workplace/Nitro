rule samtools_merge_v1_15_1:
    input:
        bam = lambda wildcards: expand(
            rules.samtools_sort_v1_15_1.output.sort, 
            sample=wildcards.sample, 
            file_id=config["SampleConfig"][wildcards.sample]
        ),
    output:
        merge = config["BaseDir"] / config["OutputDir"] / "output/{sample}/{sample}_merged.bam",
    params: 
    threads: 1
    log: 
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_merge.log"
    conda:
        "env.yaml"
    wrapper: 
        "file:rules/samtools-1.15.1/merge_wrapper.py"
