from uuid import uuid4
from pathlib import Path
from sys import argv, exit
from .banner import ShowBanner
from snakemake import snakemake
from .server import StartLogger
from rich.console import Console
from .__init__ import __version__
from argparse import ArgumentParser
from pandas import DataFrame, read_csv
from contextlib import redirect_stderr, redirect_stdout
from .qc import (
    GetWorkflowSummary,
    CheckDuplicateHashes,
    CreateOutputDirectory,
    CheckSampleSheetStructure,
    VerifyIntegrityandFileExistence,
)

# Might be removed
import yaml
import requests
from time import sleep
from threading import Thread

console = Console(tab_size=2)


def ShowArguments():
    ShowBanner()
    parser = ArgumentParser(description="Process some integers", usage="")
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=__version__,
        help="Check the version of this package",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        type=Path,
        metavar="\b",
        required=True,
        help="Provide the output directory where the results will be stored",
    )
    parser.add_argument(
        "--samplesheet",
        "-ss",
        type=Path,
        metavar="\b",
        required=True,
        help="Provide the location of samplesheet file",
    )
    # parser.add_argument("--config", "-c", metavar="\b", help="Provide the config file")
    # parser.add_argument(
    #     "--config-template",
    #     "-ct",
    #     metavar="\b",
    #     help="Get the config template that is needed by the tool",
    # )
    # parser.add_argument(
    #     "--conda-env",
    #     "-cv",
    #     metavar="\b",
    #     help="Provide folder where the conda environments will be stored",
    # )
    # parser.add_argument("--silent", "-s", metavar="\b", help="Silent the verbose")
    # parser.add_argument(
    #     "--cores", "-@", metavar="\b", help="Provide the number of cores"
    # )
    # parser.add_argument(
    #     "--breakpoints",
    #     "-bp",
    #     metavar="\b",
    #     help="Number of breakpoints required to be considered a recombinant",
    # )
    # parser.add_argument(
    #     "--parents",
    #     "-p",
    #     metavar="\b",
    #     help="Number of parents required for the sequence to be considered recombinant",
    # )
    # parser.add_argument(
    #     "--enable-deletions", "-ed", metavar="\b", help="Enable deletions"
    # )
    # parser.add_argument(
    #     "--clades", "-cl", metavar="\b", help="Provides the clade labels"
    # )
    # parser.add_argument("--update", metavar="\b", help="Update the tool")
    # parser.add_argument("--update-data", metavar="\b", help="Update the data in tool")
    # parser.add_argument(
    #     "--verify-installation", metavar="\b", help="Verify the installation"
    # )

    args = parser.parse_args()
    if len(argv) == 1:
        parser.print_help()
        exit()

    status = console.status("Checking samplesheet format")
    status.start()

    # Check the structure of samplesheet is as required
    CheckSampleSheetStructure(args.samplesheet, console, status)

    # Check if the output directory exists if not create it
    CreateOutputDirectory(args.outdir)

    # Reading the samplesheet
    status.update("Reading samplesheet")
    samplesheet = read_csv(args.samplesheet, encoding="utf-8", low_memory=False)

    # Creating configuration file
    config = {}
    config["JobID"] = uuid4().__str__()
    config["OutputDir"] = str(args.outdir)

    # Finding duplicate md5-hash
    CheckDuplicateHashes(samplesheet)

    # Rename the files and store in the output directory in the folder source
    config["SampleDict"] = VerifyIntegrityandFileExistence(
        samplesheet, console, status, args.outdir
    )

    status.update("Creating configuration")
    status.stop()

    # Getting the workflow status
    workflow, total = GetWorkflowSummary(console)

    Thread(target=StartLogger, args=[console, workflow, total]).start()
    with open("stderr.txt", "w") as f:
        with redirect_stderr(f):
            success = snakemake(
                "nitro/workflow/Snakefile",
                cores=4,
                use_conda=True,
                stats="statistics.txt",
                workdir="nitro/workflow",
                conda_prefix="condaenv",
            )
    if not success:
        try:
            requests.get("http://localhost:5000/exit")
        except:
            console.print("Connection refused because already deleted")
