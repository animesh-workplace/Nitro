rule qualimap_bamqc_v2_2_2:
    input:
        bam = rules.samtools_merge_v1_15_1.output.merge,
        index = rules.samtools_index_v1_15_1.output.index,
    output:
        directory(config["BaseDir"] / config["OutputDir"] / "output/{sample}/qualimap_report/"),
    params: 
    threads: 1
    log: 
        config["BaseDir"] / config["OutputDir"] / "output/{sample}/log/{sample}_qualimap.log"
    conda:
        "env.yaml"
    wrapper:
        "file:rules/qualimap-2.2.2/bamqc_wrapper.py"
