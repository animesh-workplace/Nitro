import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Qualimap BamQC",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
shell(
    """
        qualimap bamqc -bam {snakemake.input.bam:q} \
        -nt {snakemake.threads} -c -outdir {snakemake.output[0]:q} {log}
    """
)
requests.post(
    "http://localhost:5000/print",
    data={
        "task": "Qualimap BamQC",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}",
    },
)
