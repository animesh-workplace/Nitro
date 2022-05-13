def random(sample):
    if(sample == 'Sample1'):
        return ['1ef9b908f472643b30e3db6ba1018246']
    elif(sample == 'Sample2'):
        return [
            'd83f0ac38126c894c8fe89cebd6f4b4b',
            '9da58faf111da7918a8198845beaeafb',
        ]

rule samtools_merge_v1_15_1:
    input:
        bam = lambda wildcards: expand(
            rules.samtools_sort_v1_15_1.output.sort, 
            sample=wildcards.sample, 
            file_id=random(wildcards.sample)
        ),
    output:
        merge = "output/{sample}/{sample}_merged.bam",
    params: 
    threads: 1
    log: 
        "output/{sample}/log/{sample}_merge.log"
    conda:
        "env.yaml"
    wrapper: 
        "file:rules/samtools-1.15.1/merge_wrapper.py"
