import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Bedtools Genomecov",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
shell(
    """
        bedtools genomecov -ibam {snakemake.input.bam:q} -bga | awk -v cov=1 '$4<cov' \
        | bedtools merge -i - > {snakemake.output.low_vcf:q} {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Bedtools Genomecov",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
