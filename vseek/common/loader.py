import json
from pathlib import Path
import pandas as pd

# vseek imports
import vseek.common.vseek_paths as vsp
from vseek.common.errors import InvalidFileError
import vseek.apis.string_db as string_db
from vseek.common.io_files import *


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


def load_genes_metadata(file_path: str) -> dict:
    """returns viral genes metadata as dictioanry

    Parameters
    ----------
    file_path : str
        file path pointing to genes metadata json file

    Returns
    -------
    dict
        genes meta data
    """
    if not Path(file_path).is_file:
        raise FileNotFoundError("viral genes meta data not found")
    with open(file_path, "r") as gene_file:
        gene_metadata = json.load(gene_file)

    return gene_metadata


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
    load_path = Path(vsp.db_path()) / "final_filtered_bat_virus.csv.gz"
    if not load_path.is_file():
        raise FileNotFoundError("Unable to find viral bat database")

    return pd.read_csv(load_path)  # .drop("Unnamed: 0", axis="columns")


def load_geolocations() -> pd.DataFrame:
    """Loads country geo-locations as a pandas DataFrame

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


def load_iban_iso_codes() -> pd.DataFrame:
    """Loads alpha iso 3 codes from IBAN Database

    Returns
    -------
    pd.DataFrame
        country and their alpha ios3 codes
    """
    load_path = Path(vsp.db_path()) / "iban_iso_codes.csv"
    if not load_path.is_file():
        print("Warning: ISO files does not exists. Downloading...")
        df = pd.read_html("https://www.iban.com/country-codes")
        if len(df) == 1:
            df = df[0]
        df = df[["Country", "Alpha-2 code", "Alpha-3 code"]]
        df.columns = ["country", "iso_alpha", "alpha_iso3"]
        df.to_csv(load_path, index=False)
        return df
    else:
        return pd.read_csv(load_path)


def load_dbat_vir_db() -> pd.DataFrame:
    load_path = Path(vsp.db_path()) / "DBatVir_db.csv.gz"
    if not load_path.is_file():
        raise FileNotFoundError("Unable to find geolocations data")

    return pd.read_csv(load_path)


def load_viral_genes(accession: str) -> list[str]:
    """Returns annotated viral gene sequences

    accession : str
        accession id

    Returns
    -------
    list[str]
        list of coding sequences
    """
    viral_genome_paths = get_viral_genome_fasta_paths(query=accession)
    sel_viral_genome_path = viral_genome_paths[accession]
    header, viral_genome = load_genome(sel_viral_genome_path)

    viral_genes_paths = get_genome_genes_paths(query=accession)
    viral_genes_paths = viral_genes_paths[accession]
    meta_data = load_genes_metadata(viral_genes_paths)

    sequences = []
    for id, gene_metadata in meta_data[accession].items():

        # some genes are not annotated
        try:
            beg, end = tuple(gene_metadata["annotation"])
        except KeyError:
            continue

        annotated_sequence = viral_genome[beg : end + 1]
        sequences.append(annotated_sequence)

    return sequences


def load_species_atlas() -> pd.DataFrame:
    """Loads StringDB species atlas

    Returns
    -------
    pd.DataFrame
        taxon ids associated with species

    Raises
    ------

    """
    load_path = Path(vsp.ppi_db_path()) / "stringdb_species_codes.tsv"
    if not load_path.is_file():
        print("Warning: species atlas not found. Downlaoding...")
        species_df = string_db.download_species_atlas()
        return species_df
    else:
        return pd.read_table(load_path)


def load_human_ppi() -> pd.DataFrame:
    ppi_path = Path(vsp.ppi_db_path()) / "human_ppi.tsv.gz"
    if not ppi_path.is_file():
        print("Warning: Protein interactions data not found. Downloading ...")
        ppi_df = string_db.download_all_human_interactions()
        return ppi_df
    else:
        return pd.read_table(ppi_path, sep="\t")


# -----------------------------
# private functions (Format usage only)
# -----------------------------
def _read_fasta(fpath: str) -> str:
    """Reads FASTA file and returns it's contents

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
        # check_fasta_format(contents)

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
