# sra_callers.py
# module contains wrapper

from typing import Union
import subprocess
from pathlib import Path

# program imports
import vseek.common.vseek_paths as vsp
from vseek.common.io_files import get_prefetch_files
from vseek.common.errors import ExecutionError, SequenceFormatNotSupported
from vseek.common.checks import dependency_check


def download_fasta(
    sra_ids: Union[str, list[str]],
    threads=4,
    seq_format="fasta",
    overwrite=False,
) -> str:
    """Downloads fasta file to local machine with given sra ascension id.
    Contains 3 processes. Prefetching the data, downloading the required files
    from the NCBI's database in preparation for downloading the the fasta files.
    vbd_view step, checks for data corruption. Fasterq-dump, uses the prefetched
    files to easily download the requested fasterq-files associated with the sra
    id.

    Parameters
    ----------
    sra_ids : Union[str, list[str]]
        single string that are delimited by white space if multiple sra ids
        by white spaces or a list of sra ascension ids
    threads : int, optional
        number of threads to use for downloading fastq files,
        by default 4
    seq_format : str, optional
        Type of sequence format to be download. by default "fasta"
    overwrite : bool, optional
        Sequence files will be overwritten if set to True.
        by default False

    Returns
    ------
    str
        returns the absolute path where the fastq files are downloaded
    """

    # type checking
    if isinstance(sra_ids, str):
        sra_ids = sra_ids.split()

    # creating a prefetch directory
    prefetch_dir = vsp.init_prefetch_dir()
    fasta_dir = vsp.init_fasta_dir()

    # prefetching sra files
    srr_string = " ".join(sra_ids)
    print(f"Prefetching {srr_string} data...")
    prefetched_files = prefetch_srr_files(sra_ids, prefetch_dir)

    # validate prefetch files
    validate_prefetch_files(prefetched_files)

    # download sequence files
    print(f"Downloading {seq_format} files")
    fasterq_dump(
        prefetched_files,
        outdir=fasta_dir,
        threads=threads,
        seq_format=seq_format,
        overwrite=overwrite,
    )

    return fasta_dir


def prefetch_srr_files(accession_ids: Union[str, list[str]], out_dir: str) -> list[str]:
    """Downloads prefetched files using the prefetch executable

    Parameters
    ----------
    accession_ids : Union[str,list[str]]
        string (white spaced if multiple) or list of accession id numbers
    out_dir : str
        path to save prefetched files

    Raises
    ------
    ExecutionError
        Failed if the prefetch executable returns a error code.
    """

    # checking dependencies
    prefetch_prog = "prefetch"
    dependency_check(prefetch_prog)

    # checking types
    if isinstance(accession_ids, str):
        accession_ids = accession_ids.split()

    # prefetch cmd
    accession_ids_str = " ".join(accession_ids)
    prefetch_cmd = f"{prefetch_prog} {accession_ids_str} -O {out_dir}"
    _call(prefetch_cmd)

    # collecting all prefetch file paths
    prefetched_files = get_prefetch_files()

    return prefetched_files


def validate_prefetch_files(prefetch_file: Union[str, list[str]]) -> None:
    """Checks integrity of prefetch files using vdb-validate executable

    Parameters
    ----------
    prefetch_file : Union[str, list[str]]
        string (white spaced if multiple) or list of paths pointing to prefetch files

    Raises
    ------
    ExecutionError
        Failed if the vdb-validate executable returns a error exit code.
    """

    # checking dependencies
    validate_prog = "vdb-validate"
    dependency_check(validate_prog)

    # type checking
    if isinstance(prefetch_file, str):
        prefetch_file = prefetch_file.split()

    for p_file_path in prefetch_file:
        cmd = f"{validate_prog} {p_file_path}"
        _call(cmd)


def fasterq_dump(
    prefetch_file: Union[str, list[str]],
    outdir: str,
    seq_format="fasta",
    threads=4,
    overwrite=False,
) -> None:
    """Downloads sequence files ("Fasta" or "Fastq") using the fasterq-dump executable

    Parameters
    ----------
    prefetch_file : Union[str, list[str]]
        string (white spaced if multiple) or list of prefetch file paths.
    outdir : str
        directory to store downloaded fasta files
    threads : int, optional
        number of threads used to download files, by default 4
    seq_format : str
        file format output. supported formats = ["fastq", "fasta"]
    overwrite : bool, optional
        Sequence files will be re-downloaded and overwritten if set to True.
        by default False

    Raises
    ------
    ExecutionError
        _description_
    """

    # checking dependencies
    fastq_prog = "fastq-dump"
    fasterq_prog = "fasterq-dump"
    dependency_check(fastq_prog)
    dependency_check(fasterq_prog)

    # type checking
    if isinstance(prefetch_file, str):
        prefetch_file = prefetch_file.split()

    # downloading sequence files
    for prefetch_path in prefetch_file:

        # checking if the file exists. If overwrite
        name = prefetch_path.rsplit("/", 1)[-1].split(".")[0]
        if overwrite is False:
            path_obj = Path(vsp.metagenome_path()) / f"{name}.fasta"
            check = path_obj.is_file()
            if check is True:
                print(f"{name}.fasta already exists... skipping")
                continue

        if seq_format == "fastq":
            fasterq_cmd = f"{fasterq_prog} {prefetch_path} -e {threads} -O {outdir}"
            _call(fasterq_cmd)
        elif seq_format == "fasta":
            fastq_cmd = f"{fastq_prog} {prefetch_path} --{seq_format} 60 -O {outdir}"
            _call(fastq_cmd)
        else:
            raise SequenceFormatNotSupported(f"{seq_format} is not supported")


# -----------------------------
# Private functions
# -----------------------------
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
        Failed if the executable returns a error exit code.
    """
    cmd = cmd.split()
    call = subprocess.run(
        cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if call.returncode != 0:
        # TODO: add a parser here to print the error since it is captured in stderr
        raise ExecutionError("Execution failed. Non-zero exit status captured")
