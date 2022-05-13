import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Samtools Sort",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
shell(
    """
        samtools sort -@ {snakemake.threads} {snakemake.input[0]:q} -o {snakemake.output.sort:q} {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Samtools Sort",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
