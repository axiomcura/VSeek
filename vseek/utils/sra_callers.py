import shutil
import subprocess
from pathlib import Path

# program imports
from vseek.common.errors import ExecutionError
from vseek.common.checks import dependency_check


def download_fastq(
    sra_ids: str | list,
    threads=4,
    prefetch_dir="SRA_download",
    fastq_dir="fastq_files",
    verbose=False,
) -> str:
    """Downloads fastq file to local machine with given sra ascension id.
    Contains 3 processes. Prefetching the data, downloading the required files
    from the NCBI's database in preperation for downloading the the fastq files.
    vbd_view step, checks for data corruption. Fasterq-dump, uses the prefetched
    files to easily download the requested fasterq-files associated with the

    Parameters
    ----------
    sra_ids : str | list
        single string that are delimited if multiple sra ids
        by white spaces or a list of sra ascension ids
    threads : int, optional
        number of threads to use for downloading fastq files,
        by default 4
    prefetch_dir : str, optional
        directory for storing prefetched sra files,
        by default "SRA_download"
    fastq_dir : str, optional
        directory name for saving fastq files
        by default "fastq_files"

    Returns
    ------
    str
        returns the absolute path where the fastq files are downloaded
    """

    # type checking
    if isinstance(sra_ids, str):
        sra_ids = sra_ids.split()

    # creating a prefetch directory
    prefetch_path = Path(f"results/{prefetch_dir}")
    prefetch_path.mkdir(exist_ok=True)
    fastq_path = Path(f"results/{fastq_dir}")
    fastq_path.mkdir(exist_ok=True)

    # executeable names
    prefetch_prog = "prefetch"
    fasterq_prog = "fasterq-dump"

    # checking dependencies
    dependency_check(prefetch_prog)
    dependency_check(fasterq_prog)

    # prefetching sra files
    print("Prefetching SRR data...")
    sra_ids_str = " ".join(sra_ids)
    prefetch_cmd = f"{prefetch_prog} {sra_ids_str} -O {prefetch_path.absolute()}"
    _call(prefetch_cmd)

    # TODO: add vdb_view steps checking that the files are not corrupt
    # downloading fastq files
    print("Downloading Fastq files")
    for sra_id in sra_ids:
        fasterq_cmd = f"{fasterq_prog} {prefetch_path.absolute()}/{sra_id}/{sra_id}.sra -e {threads} -O {fastq_path.absolute()}"
        print(fasterq_cmd)
    #    _call(fasterq_cmd)
    #    print(f"{sra_id} fastq file download complete")

    return fastq_path.absolute()


# TODO: add error parser for each callable?
def _call(cmd: str) -> int:
    """Wrapper for calling executable

    Parameters
    ----------
    cmd : str
        command input

    Returns
    -------
    int
        Return code


    Raises:
    -------
    ExecutionFailedError:
        raised when a non-zero exit status code is raised from the executable
    """
    cmd = cmd.split()
    call = subprocess.run(
        cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if call.returncode != 0:
        # TODO: add a parser here to print the error since it is captured in stderr
        raise ExecutionError("Execution failed. Non-zero exit status captured")
