import os
from pathlib import Path

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
    db_path = str(root / "db")

    return db_path


def genome_db_path() -> str:
    """Returns database path

    Returns
    -------
    str
        path to genome database
    """
    gdb_path = Path(db_path()) / "genome"
    return str(gdb_path.absolute())


# -----------------------------
# Initializations of databases
# -----------------------------
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
