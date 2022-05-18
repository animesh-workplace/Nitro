import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "BCFTools Consensus",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
shell(
    """
        bcftools consensus -f {snakemake.params.reference:q} {snakemake.input.vcf:q} -m {snakemake.input.mask:q} -o {snakemake.output[0]:q} {log}
    """
)
# bcftools maskfasta -fi {snakemake.output[0]:q} -fo Sample1.consensus_masked.fasta -bed {snakemake.input.no_vcf:q} {log}
# sed -i "1s/.*/>{snakemake.wildcards[0]}/" {snakemake.output[0]:q}
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "BCFTools Consensus",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
