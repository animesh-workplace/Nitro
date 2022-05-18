import re
import json
import shutil
from hashlib import md5
from pathlib import Path
from pandas import read_csv
from itertools import groupby
from collections import Counter
from snakemake import snakemake
from snakemake.shell import shell
from mmap import mmap, ACCESS_READ
from importlib.util import find_spec
from contextlib import redirect_stdout
from tempfile import TemporaryDirectory

# from pydot import graph_from_dot_file


def CheckSampleSheetStructure(samplesheet, console, status):
    # Check the structure of samplesheet is as required
    tempdir = TemporaryDirectory()
    BASE_DIR = Path(find_spec("nitro").origin).parent
    try:
        shell(
            """
                chkcsv.py -r {samplesheet:q} \
                -f {BASE_DIR:q}/formatspec/samplesheet.fmt 2> {tempdir.name:q}/csv.log
            """
        )
    except:
        status.stop()
        console.rule("Error", style="red")
        with open(f"{tempdir.name}/csv.log") as f:
            console.print(f.read().replace("Error: ", ""), style="red", highlight=False)
        console.rule(style="red")

    console.print("[green]:heavy_check_mark: Verified the samplesheet")


def CreateOutputDirectory(outdir):
    """
    Check if the output directory exists,
    if yes delete it and then recreate the folder else
    if create it
    """
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
        # This source directory holds renamed fastq file for workflow
        p = outdir / "source"
        p.mkdir(parents=True, exist_ok=True)


def CheckDuplicateHashes(samplesheet):
    # Finding duplicate md5-hash
    duplicate_files = samplesheet[samplesheet.duplicated(["MD5-hash"], keep=False)].File
    if len(duplicate_files):
        console.print("Duplicate hash found for IDS:", style="red", highlight=False)
        for file in duplicate_files:
            console.print(f"\t- {file}", style="red", highlight=False)
        status.stop()
        exit()


def VerifyIntegrityandFileExistence(samplesheet, console, status, outdir):
    # Rename the files and store in the output directory in the folder source
    status.update(f"Verifing file integrity (0/{len(samplesheet)})")
    hash_failed = 0
    exist_failed = []
    config = {}
    for index, row in samplesheet.iterrows():
        InputFile = Path(row.File)
        if row.Sample not in config:
            config[row.Sample] = []
        config[row.Sample].append(row["MD5-hash"])
        # If InputFile doesnot exist then print and continue checking the next
        if not InputFile.exists():
            exist_failed.append(f"{InputFile} (Not found)")
            continue

        with open(InputFile) as file, mmap(
            file.fileno(), 0, access=ACCESS_READ
        ) as file:
            calculated_hash = md5(file).hexdigest()
            if calculated_hash == row["MD5-hash"]:
                status.update(f"Verifing file integrity ({index+1}/{len(samplesheet)})")
            else:
                console.print(
                    f"File integrity failed for {row.File}",
                    style="red",
                    highlight=False,
                )
                hash_failed = hash_failed + 1
        suffix = "".join(InputFile.suffixes)

        """
            Checking whether the file is gunzipped or not,
            if not then gunzip the file
        """
        FilePath = outdir / "source" / f"{row.Sample}_{calculated_hash}{suffix}"

        if InputFile.suffix == ".gz":
            if not FilePath.exists():
                shell(
                    """
                        cp {row.File} {FilePath:q}
                    """
                )
        else:
            FilePathgz = (
                outdir / "source" / f"{row.Sample}_{calculated_hash}{suffix}.gz"
            )
            if not FilePathgz.exists():
                shell(
                    """
                        cp {row.File} {FilePath:q}
                        gzip {FilePath:q}
                    """
                )

    # Guard check for file existence failure
    if exist_failed:
        console.print("[red]✖ Failed file existence check")
        for file in exist_failed:
            console.print(f"\t- {file}", style="red", highlight=False)
        status.stop()
        exit()

    # Guard check for file integrity failure
    if hash_failed:
        console.print("[red]✖ Failed file integrity")
        status.stop()
        exit()

    console.print(
        f"[green]:heavy_check_mark: Verified the file hash ({len(samplesheet)}/{len(samplesheet)})"
    )
    return config


def GetWorkflowSummary(console, config_loc):
    tempdir = TemporaryDirectory()
    summary = Path(tempdir.name, "summary.json")

    # Rename rule schema
    rule_schema = {
        "all": "All",
        "fastp": "FastP",
        "fastqc": "FastQC",
        "bwa_mem": "BWA MEM",
        "samtools_sort": "Samtools Sort",
        "trimmomatic_se": "Trimmomatic SE",
        "qualimap_bamqc": "Bedtools Genomecov",
        "bedtools_genomecov": "Qualimap BamQC",
        "bedtools_maskfasta": "Bedtools Maskfasta",
        "samtools_merge": "Samtools Merge",
        "samtools_index": "Samtools Index",
        "gatk_haplotype": "GATK HaplotypeCaller",
        "bcftools_consensus": "BCFTools Consensus",
    }

    """
        Getting the summary stats in a file in temporary directory
        using the dry run of snakemake that emits a table
    """
    with open(summary, "w") as f:
        with redirect_stdout(f):
            snakemake(
                "nitro/workflow/Snakefile",
                cores=4,
                quiet=True,
                dryrun=True,
                nocolor=True,
                use_conda=True,
                conda_prefix="condaenv",
                configfiles=[config_loc],
                workdir="nitro/workflow",
            )

    with open(summary, "r") as f:
        sample = f.read()

    all_input_lines = sample.splitlines()
    """
        Guard clause if the snakemake dry run is showing no output,
        meaning that the output files are already present
    """
    if len(all_input_lines) == 0:
        console.print(
            "\t[italic][yellow]Nothing to be done (all requested files are present and up to date)",
            highlight=False,
        )
        exit()
    all_input_lines = all_input_lines[1:]

    # Converting the table into dictionary
    row_dict = []
    for is_blank, input_lines_iter in groupby(
        all_input_lines, key=lambda s: not bool(s.strip())
    ):
        input_lines = list(input_lines_iter)
        if is_blank:
            continue

        # First two lines are field names and separator dashes
        names, dashes = input_lines[:2]

        # Using regex to get start/end of each '---' divider, and making slices
        spans = (match.span() for match in re.finditer("-+", dashes))
        slices = [slice(sp[0], sp[1] + 1) for sp in spans]
        # Slice and get the row name
        names = [names[sl].strip() for sl in slices]

        for line in input_lines[2 : len(input_lines) - 1]:
            temp = dict(zip(names, (line[sl].strip() for sl in slices)))
            temp["job"] = rule_schema[temp["job"].split("_v")[0].replace("_main", "")]
            if temp["job"] != "All":
                row_dict.append(temp)

    # Counting the total jobs
    total = sum([int(row["count"]) for row in row_dict])

    # Creating a unified dictionary for helping in printing
    job_count = {
        row["job"]: {"total": int(row["count"]), "running": 0, "completed": 0}
        for row in row_dict
    }

    return job_count, total
