import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "FastP",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
shell(
    """
        fastp --thread {snakemake.threads} -i {snakemake.input[0]:q} \
        {snakemake.params} -o {snakemake.output.trim:q} \
        --json {snakemake.output.json:q} --html {snakemake.output.html:q} {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "FastP",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
