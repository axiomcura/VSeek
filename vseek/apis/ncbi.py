import sys
from time import sleep
from io import StringIO
from typing import Union
from pathlib import Path

import requests
import pandas as pd
import entrezpy.conduit

# VSeek imports
import vseek.common.vseek_paths as vsp
from vseek.common.checks import genome_db_exist
from vseek.utils.parsers import parse_ncbi_viral_accessions
from vseek.common.io_files import genome_db_files, save_genome


def get_all_viral_accessions() -> pd.DataFrame:
    """Sends a request to ncbi to obtain all viral genome accessions"""

    csv_name = "ncbi_viral_accession_db.csv.gz"
    db = vsp.init_db_path()
    file_path_obj = Path(db) / csv_name
    check = file_path_obj.is_file()
    if not check:
        print("downloading NCBI's viral genome accession database")
        url = "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            raise ConnectionError("could not donwload data")

        # parsing contents
        header, viral_cont = parse_ncbi_viral_accessions(r.text)

        # making a dataframe and saving it
        save_path = str(file_path_obj.absolute())
        ncbi_viral_df = pd.DataFrame(data=viral_cont, columns=header)

        # saving file
        ncbi_viral_df.to_csv(save_path, compression="gzip", index=False)

        return ncbi_viral_df

    else:
        file_path = str(file_path_obj.absolute())
        return pd.read_csv(file_path)


def get_viral_genomes(email: str, accessions: Union[str, list], buffer=0.5) -> dict:
    """Downloads all viral gneomes and generates a viral genome profiles

    Parameters
    ----------
    email : str
        valid email address require to send requests to the NCBI
        database using entrez.
    accession : Union[str, list]
        string or list of accession numbers
    buffer : int, float
        Buffer time added after submitting a request in seconds
        Default = 0.5

    Return
    ------
    None
        Generates a genome database under ./db
    """
    print("\nBuilding genome database ...")
    if isinstance(accessions, str):
        accessions = accessions.split()

    # if the database exists, check for missing files and download them
    print("Checking if genome database exists ...")
    if genome_db_exist():
        print("Genome database already exists. Checking for missing files")
        expected = set(accessions)
        gdb_files = {gpath.rsplit("/", 1)[-1] for gpath in genome_db_files()}
        missing = expected - gdb_files

        # download missing genomes
        if len(missing) > 0:
            print("Warning there are some missing genomes")
            for acc_id in missing:
                print(f"> Requesting {acc_id} genome ... ")
                viral_genome = _call_entrez_viral_genome(
                    email=email, accession=acc_id, buffer=buffer
                )
                save_genome(accession=acc_id, contents=viral_genome)

            return genome_db

        # no missing files, just return path
        else:
            print("No missing genomes")
            return genome_db

    # if database does not exist, create database and download all genomes
    genome_db = vsp.init_genome_db_path()
    for acc_id in accessions:
        print(f">Requesting {acc_id} genome ... ")
        viral_genome = _call_entrez_viral_genome(
            email=email, accession=acc_id, buffer=buffer
        )
        save_genome(accession=acc_id, contents=viral_genome)

    return genome_db


# ------------------------------
# private caller functions
# ------------------------------
def _call_entrez_viral_genome(email: str, accession: str, buffer=0.3) -> str:
    """Submits request to NCBI viral genome database via entrez portal.

    Parameters
    ----------
    email : str
        valid email address
    accession : str
        genome accession number
    buffer : int, float
        Buffer time added when submitting a request in seconds
        Default = 0.3

    Returns
    -------
    str
        Genome fasta
    """
    # creating a reference of the default stdout
    old_stdout = sys.stdout

    # creating a container storing stdout
    fasta_result = StringIO()
    sys.stdout = fasta_result

    # calling
    c = entrezpy.conduit.Conduit(email)
    fetch_influenza = c.new_pipeline()
    sid = fetch_influenza.add_search(
        {
            "db": "nucleotide",
            "term": f"{accession}",
            "rettype": "count",
            "datetype": "pdat",
        }
    )
    fid = fetch_influenza.add_fetch(
        {"retmax": 10, "retmode": "text", "rettype": "fasta"}, dependency=sid
    )
    c.run(fetch_influenza)
    sleep(buffer)

    # store the string from the stdout into variable
    genome = fasta_result.getvalue()

    # now redict back stdout to screen
    sys.stdout = old_stdout

    return genome
