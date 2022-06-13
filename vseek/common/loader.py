import glob
from pathlib import Path
import pandas as pd

# vseek imports
import vseek.common.vseek_paths as vsp
from vseek.common.checks import check_fasta_format
from vseek.common.errors import InvalidFileError

# TODO: implement
def load_profile() -> dict:
    """Load viral genome profile

    Returns
    -------
    dict
        viral genome profile
    """
    pass


def load_genome(file_path: str) -> tuple:
    """

    Parameters
    ----------
    file_path : str
        path to fasta file

    Returns
    -------
    tuple
        header, flat sequence in a tuple
    """

    # extension check
    contents = _read_fasta(file_path)

    # flattening sequence
    header, sequence = _flatten_fasta_sequence(contents)
    return (header, sequence)


def load_bat_virus_data() -> pd.DataFrame:
    """Returns bat virus database

    Returns
    -------
    pd.DataFrame
        bat virus database

    Raises
    ------
    FileNotFoundError
        raise if the file is not found in the database
    """
    load_path = Path(vsp.db_path()) / "filtered_bat_virus.csv.gz"
    if not load_path.is_file():
        raise FileNotFoundError("Unable to find viral bat database")

    return pd.read_csv(load_path)


def load_geolocations() -> pd.DataFrame:
    """Loads country geo-locations as a pandas dataframe

    Returns
    -------
    pd.DataFrame
        geo locations

    Raises
    ------
    FileNotFoundError
        raised file is not found in the database
    """
    load_path = Path(vsp.db_path()) / "geolocations.csv.gz"
    if not load_path.is_file():
        raise FileNotFoundError("Unable to find geolocations data")

    return pd.read_csv(load_path)


# -----------------------------
# private functions (Format usage only)
# -----------------------------
def _read_fasta(fpath: str) -> str:
    """Reads FASTA file and returns it's contens

    Parameters
    ----------
    fpath : str
        file path that points to fasta file

    Returns
    -------
    str
        Fasta file contents
    """
    if not fpath.endswith(".fasta"):
        raise InvalidFileError("Fasta files are either .fasta or .fa")

    with open(fpath, "r") as infile:
        contents = infile.read()

        # check if the contents
        check_fasta_format(contents)

    return contents


def _flatten_fasta_sequence(contents: str) -> str:
    """Flattens FASTA sequence into one single line

    Parameters
    ----------
    contents : str
        original fasta contents

    Returns
    -------
    tuple
        contains header and flatten sequence in a tuple
    """
    split_contents = contents.splitlines()
    header = split_contents[0]
    sequence = "".join(split_contents[1:])

    return (header, sequence)
