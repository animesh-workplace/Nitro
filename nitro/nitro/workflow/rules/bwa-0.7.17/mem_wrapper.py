import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "BWA MEM",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
shell(
    """
        bwa mem {snakemake.params.reference:q} {snakemake.input[0]:q} -t {snakemake.threads} \
        -R '@RG\\tID:fcl_1_lane1\\tPL:Illumina\\tSM:{snakemake.wildcards[0]}' \
        -T 1 -M -o {snakemake.output.align:q} {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "BWA MEM",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
