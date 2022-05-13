import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "GATK HaplotypeCaller",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
shell(
    """
        gatk --java-options "-Xmx4g" HaplotypeCaller -R {snakemake.params.reference:q} \
        -I {snakemake.input.bam:q} -O {snakemake.output[0]:q} {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "GATK HaplotypeCaller",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
