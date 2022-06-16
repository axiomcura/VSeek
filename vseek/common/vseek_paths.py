import os
from pathlib import Path

# vseek imports
from vseek.common.checks import results_dir_exists

# -----------------------------
# Path functions
# -----------------------------
def relative_root_path() -> str:
    re_path = "./VSeek"
    return re_path


def root_path() -> str:
    """Returns project root path

    Returns
    -------
    str
        project root path
    """
    splitter = "VSeek"
    working_dir = os.getcwd()
    root_path = working_dir.split(splitter)[0] + splitter

    return root_path


def db_path() -> str:
    """Returns path pointing to database folder

    Returns
    -------
    str
        path to db directory
    """
    root = Path(root_path())
    database_path = str((root / "db").absolute())

    return database_path


def results_dir() -> str:
    """Returns a path pointing to

    Returns
    -------
    str
        _description_
    """
    root = Path(root_path())
    results_path = str((root / "results").absolute())

    return results_path


def genome_db_path() -> str:
    """Returns database path

    Returns
    -------
    str
        path to genome database
    """
    gdb_path = Path(db_path()) / "genome"
    return str(gdb_path.absolute())


def prefetch_path() -> str:
    """Returns path to SRA Prefetch directory

    Returns
    -------
    str
        path to SRA Prefetch directory
    """
    # checking if the results directory exists if not create one
    if not results_dir_exists():
        init_results_dir()

    pfetch_path = Path(results_dir()) / "SRA_prefetch"
    return str(pfetch_path.absolute())


def metagenome_path() -> str:
    """Returns path to fasta meta-genome directory

    Returns
    -------
    str
        path to fasta directory
    """
    # checking if the results directory exists if not create one
    if not results_dir_exists():
        init_results_dir()

    fasta_dir_path = Path(results_dir()) / "fasta_files"
    return str(fasta_dir_path.absolute())


def ppi_db_path() -> str:
    """Returns path to string_db directory

    Returns
    -------
    str
        path to string db directory
    """
    ppi_path = Path(db_path()) / "protein_interactions"
    return str(ppi_path.absolute())


# -----------------------------
# Initializations of directories
# -----------------------------
def init_results_dir() -> str:
    """Creates results directory if it does not exists.

    Returns
    -------
    str
        path to results dictionary
    """
    results_path_obj = Path(results_dir())
    results_path_obj.mkdir(exist_ok=True)

    return str(results_path_obj.absolute())


def init_db_path() -> str:
    """Creates a database folder in the root project directory.

    Returns
    -------
    str
        _description_
    """
    db_path_obj = Path(db_path())
    db_path_obj.mkdir(exist_ok=True)

    return str(db_path_obj.absolute())


def init_genome_db_path() -> str:
    """Returns path to viral genome directory. If the directory does not exists,
    it will create one.

    Returns
    -------
    str
        path to genome data
    """

    genome_path_obj = Path(genome_db_path())
    genome_path_obj.mkdir(exist_ok=True)
    genome_path = str(genome_path_obj.absolute())

    return genome_path


def init_prefetch_dir() -> str:
    """Creates a prefetch directory if it does not exists and returns
    the path.

    Returns
    -------
    str
        path to prefetch directory
    """
    pfetch_path_obj = Path(prefetch_path())
    pfetch_path_obj.mkdir(exist_ok=True)
    pfetch_path = str(pfetch_path_obj.absolute())

    return pfetch_path


def init_fasta_dir() -> str:
    """Creates a fasta directory if it does not exists and returns
    the path.

    Returns
    -------
    str
        path to fasta directory
    """
    metagenome_path_obj = Path(metagenome_path())
    metagenome_path_obj.mkdir(exist_ok=True)
    fasta_path_str = str(metagenome_path_obj.absolute())

    return fasta_path_str


def init_string_dir() -> str:
    """Creates a fasta directory if it does not exists and returns
    the path.

    Returns
    -------
    str
        path to fasta directory
    """
    ppi_path_obj = Path(ppi_db_path())
    ppi_path_obj.mkdir(exist_ok=True)
    fasta_path_str = str(ppi_path_obj.absolute())

    return fasta_path_str


def init_profile_dir() -> str:
    """Creates a profile directory if it does not exists and returns
    the path.

    Returns
    -------
    str
        path to profile directory
    """
    profile_path = Path(results_dir()) / "ppi-profiles"
    profile_path.mkdir(exist_ok=True)

    profile_path_str = str(profile_path.absolute())

    return profile_path_str


def init_plots_dir() -> str:
    """Creates a plots directory in the results directory

    Returns
    -------
    str
        path to plots directory
    """
    plots_path_obj = Path(results_dir()) / "plots"
    plots_path_obj.mkdir(exist_ok=True)

    plots_path_str = str(plots_path_obj.absolute())
    return plots_path_str
