import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Samtools Merge",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
shell(
    """
        samtools merge {snakemake.input.bam:q} -o {snakemake.output.merge:q} -@ {snakemake.threads} {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Samtools Merge",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
