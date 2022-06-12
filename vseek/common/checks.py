from pathlib import Path
import shutil

# program imports
import vseek.common.vseek_paths as vsp
from vseek.common.errors import *

# dowloading the data module
def dependency_check(prog: str) -> bool:
    """Checks if the executable program is available

    Parameters
    ----------
    prog : str
        name of the program or executable

    Returns
    -------
    Bool
        Returns True if the program exists

    Raises:
    -------
    MissingDependencyError:
        Raised If the program is not installed to local machine
    UnsupportedDependencyError:
        Raised if attempt to call unsupported executable. Please look at supported
        executables in the README.md file.
    """
    callables = ["prefetch", "fastq-dump", "fasterq-dump", "vdb-validate"]

    # checking if the callable is supported
    if prog not in callables:
        raise UnsupportedDependencyError(
            f"{prog} is not supported. Supported programs {callables}"
        )

    # checking if supported callable is installed
    check = shutil.which(prog)
    if check is None:
        raise MissingDependencyError(f"{prog} is not installed")

    return True


def check_fasta_format(fasta_content: str) -> None:
    """Checks the integrity of the FASTA format

    Parameters
    ----------
    fasta_content : str
        Fasta file content

    Raises
    ------
    InvalidFastaFormatError
        Raised if the contents possess invalid FASTA structure
    """
    fasta_content = fasta_content.splitlines()
    header = fasta_content[0]

    # header check
    if not header.startswith(">"):
        raise InvalidFastaFormatError("Loaded fasta file contains a corrupt header")

    # sequence length (each line should be 60, last line can be lower than 60)
    for line in fasta_content[1:]:
        if not len(line) <= 60:
            raise InvalidFastaFormatError(
                "Loaded fasta file contains invalid sequence length"
            )

    return


# ------------------------------
# database checks
# ------------------------------
def genome_db_exist() -> bool:
    """Checks if the genome database exists

    Returns
    -------
    bool
        True if exists, False it does not exist
    """

    genome_db_path = vsp.genome_db_path()
    if not Path(genome_db_path).is_dir:
        return False
    return True


def prefetch_dir_exists() -> bool:
    """Checks if the prefetch directory exists

    Returns
    -------
    bool
        True if exists, False it does not exist
    """
    prefetch_dir = vsp.prefetch_path()
    if not Path(prefetch_dir).is_dir:
        return False
    return True


def results_dir_exists() -> bool:
    """Checks if the results directory exists

    Returns
    -------
    bool
        True if exists, False it does not exist
    """
    results_dir = vsp.results_dir()
    if not Path(results_dir).is_dir():
        return False
    return True
