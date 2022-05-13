import requests
from os import path
from datetime import datetime
from snakemake.shell import shell
from tempfile import TemporaryDirectory

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "FastQC",
        "status": "Started",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)

tempdir = TemporaryDirectory()
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
zip_path = path.join(
    tempdir.name,
    f"{snakemake.wildcards.sample}_{snakemake.wildcards.file_id}_fastqc.zip",
)
html_path = path.join(
    tempdir.name,
    f"{snakemake.wildcards.sample}_{snakemake.wildcards.file_id}_fastqc.html",
)
shell(
    """
        fastqc {snakemake.params} -t {snakemake.threads} \
        --outdir {tempdir.name:q} {snakemake.input[0]:q} {log}
    """
)
if snakemake.output.html != html_path:
    shell("mv {html_path:q} {snakemake.output.html:q}")
if snakemake.output.zip != zip_path:
    shell("mv {zip_path:q} {snakemake.output.zip:q}")

requests.post(
    "http://localhost:5000/print",
    data={
        "task": "FastQC",
        "status": "Finished",
        "time": str(datetime.now()),
        "db_loc": snakemake.config["db_loc"],
        "sample": f"{snakemake.wildcards.sample}-{snakemake.wildcards.file_id}",
    },
)
