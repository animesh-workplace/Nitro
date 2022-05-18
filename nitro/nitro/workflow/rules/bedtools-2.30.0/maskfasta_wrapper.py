import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Bedtools Maskfasta",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
shell(
    """
        bedtools maskfasta -fi {snakemake.params.reference:q} -fo {snakemake.output[0]:q} -bed {snakemake.input.no_vcf:q} {log}
        sed -i "1s/.*/>{snakemake.wildcards[0]}/" {snakemake.output[0]:q}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Bedtools Maskfasta",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
