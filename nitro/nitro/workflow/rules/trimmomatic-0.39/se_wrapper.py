import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Trimmomatic SE",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
shell(
    """
        trimmomatic SE -threads {snakemake.threads} -phred33 {snakemake.input[0]:q} {snakemake.output.trim:q} \
        ILLUMINACLIP:{snakemake.params.adapters:q}:2:30:30 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:30 {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Trimmomatic SE",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
