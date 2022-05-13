import json
import shutil
from hashlib import md5
from pathlib import Path
from pandas import read_csv
from collections import Counter
from snakemake import snakemake
from snakemake.shell import shell
from mmap import mmap, ACCESS_READ
from importlib.util import find_spec
from contextlib import redirect_stdout
from tempfile import TemporaryDirectory


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
    if outdir.exists():
        shutil.rmtree(outdir)
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
    config = {"SampleName": [], "FileID": []}
    for index, row in samplesheet.iterrows():
        InputFile = Path(row.File)
        config["SampleName"].append(row.Sample)
        config["FileID"].append(row["MD5-hash"])

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
        if InputFile.suffix == ".gz":
            shell(
                """
                    cp {row.File} {outdir:q}/source/{row.Sample}_{calculated_hash}{suffix}
                """
            )
        else:
            shell(
                """
                    cp {row.File} {outdir:q}/source/{row.Sample}_{calculated_hash}{suffix}
                    gzip {outdir:q}/source/{row.Sample}_{calculated_hash}{suffix}
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

    console.print("[green]:heavy_check_mark: Verified the file hash")
    return config


def GetWorkflowSummary(console):
    tempdir = TemporaryDirectory()
    summary = Path(tempdir.name, "summary.json")

    """
        Getting the summary stats in a file in temporary directory
        using the d3dag of snakemake that emits a JSON
    """
    with open(summary, "w") as f:
        with redirect_stdout(f):
            snakemake(
                "nitro/workflow/Snakefile",
                cores=4,
                dryrun=True,
                quiet=True,
                nocolor=True,
                use_conda=True,
                printd3dag=True,
                conda_prefix="condaenv",
                workdir="nitro/workflow",
                # conda_base_path="" [donot delete will be required in the future]
            )

    # Reading back the JSON file and saving it in a dict format
    with open(summary) as f:
        temp = json.loads(f.read())
    # Rename rule schema
    rule_schema = {
        "all": "All",
        "fastp": "FastP",
        "fastqc": "FastQC",
        "bwa_mem": "BWA MEM",
        "samtools_sort": "Samtools Sort",
        "trimmomatic_se": "Trimmomatic SE",
        "qualimap_bamqc": "Qualimap BamQC",
        "samtools_merge": "Samtools Merge",
        "samtools_index": "Samtools Index",
        "gatk_haplotype": "GATK HaplotypeCaller",
        "bcftools_consensus": "BCFTools Consensus",
    }

    # console.log(temp)
    ordering = list(
        set(
            [
                rule_schema[i["value"]["rule"].split("_v")[0].replace("_main", "")]
                for i in temp["nodes"]
            ]
        )
    )
    job_count = dict(
        Counter(
            [
                rule_schema[i["value"]["rule"].split("_v")[0].replace("_main", "")]
                for i in temp["nodes"]
            ]
        )
    )
    # Counting the total jobs
    total = sum([value for key, value in job_count.items()])
    # Creating a unified dictionary for helping in printing
    job_count = {
        key: {"total": value, "running": 0, "completed": 0}
        for key, value in job_count.items()
        if (not key == "All")
    }

    return job_count, total