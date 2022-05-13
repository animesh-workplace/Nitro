rule fastqc_main_v0_11_9:
    input:
        "data/{sample}_{file_id}.fastq.gz",
    output:
        zip = "output/{sample}/fastqc_report/{sample}_{file_id}_fastqc.zip",
        html = "output/{sample}/fastqc_report/{sample}_{file_id}_fastqc.html",
    params:
    version: "0.11.9"
    threads: 1
    log:
        "output/{sample}/log/{sample}_{file_id}_fastqc.log"
    conda:
        "env.yaml"
    wrapper:
        "file:rules/fastqc-0.11.9/main_wrapper.py"
