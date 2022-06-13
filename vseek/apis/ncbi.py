import sys
from time import sleep
from io import StringIO
from typing import Union
from pathlib import Path
from collections import defaultdict

import requests
import pandas as pd
import entrezpy.conduit

# VSeek imports
import vseek.common.vseek_paths as vsp
from vseek.common.checks import genome_db_exist
from vseek.utils.parsers import parse_ncbi_viral_accessions
from vseek.common.io_files import (
    genome_dir_paths,
    get_genome_genes_paths,
    save_genome,
    save_genes,
)


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
            raise ConnectionError("could not download data")

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
    genome_db = vsp.init_genome_db_path()
    print("\nBuilding genome database ...")
    if isinstance(accessions, str):
        accessions = accessions.split()

    # if the database exists, check for missing files and download them
    print("Checking if genome database exists ...")
    if genome_db_exist():
        print("Genome database already exists. Checking for missing files")
        expected = set(accessions)
        gdb_files = {gpath.rsplit("/", 1)[-1] for gpath in genome_dir_paths()}
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


def get_viral_genes(email: str, accession: Union[str, list[str]], buffer=0.5):
    """Sends requests to obtain gene information with the associated genome

    Parameters
    ----------
    email : Union[str, list[str]]
        valid email address require to send requests to the NCBI
        database using entrez.
    accession : Union[str, list]
        string or list of accession numbers
    buffer : int, float
        Buffer time added after submitting a request in seconds
        Default = 0.5
    """
    print("\nDownloading viral genes metadata")
    # type checking
    if isinstance(accession, str):
        accession = accession.split()

    # checking of the ge
    print("Checking if genome database exists ...")
    if genome_db_exist():
        print("Database exists! Checking for missing gene data files")

        # checking for missing file
        expected = set(accession)
        current_files = {
            p.rsplit("/", 1)[-1].split("_genes")[0]
            for p in get_genome_genes_paths().values()
        }
        missing = expected - current_files

        if len(missing) > 0:
            print("Warning there are some missing genomes")
            for acc_id in missing:
                print(f"> Requesting {acc_id} genes ... ")
                viral_genome_genes = _call_entrez_viral_genes(
                    email=email, accession=acc_id, buffer=buffer
                )
                save_genes(accession=acc_id, contents=viral_genome_genes)


def generate_viral_genome_profile(genome_db: str) -> dict:
    """Generates a json file that profiles

    Parameters
    ----------
    genome_db : str
        path to genome database

    Returns
    -------
    dict
        viral genome profiles. Also written in JSON format in the genome database
    """
    # gene_positions = _call_entrez_viral_genes(email=email, accession=acc_id)
    pass


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
    fetch_genome = c.new_pipeline()
    sid = fetch_genome.add_search(
        {
            "db": "nucleotide",
            "term": f"{accession}",
            "rettype": "count",
            "datetype": "pdat",
        }
    )
    fid = fetch_genome.add_fetch(
        {"retmax": 10, "retmode": "text", "rettype": "fasta"}, dependency=sid
    )
    c.run(fetch_genome)
    sleep(buffer)

    # store the string from the stdout into variable
    genome = fasta_result.getvalue()

    # now redict back stdout to screen
    sys.stdout = old_stdout

    return genome


def _call_entrez_viral_genes(email: str, accession: str, buffer=0.5) -> dict:
    """Submits request to NCBI's genes database via entrez portal

    Parameters
    ----------
    email : str
        valid email address
    accession : str
        genome accession number
    buffer : int, float
        Buffer time added when submitting a request in seconds
        Default = 0.5

    Returns
    -------
    dict
        gene metadata
    """

    # call genes data

    # internal parser that returns a tuple of ranges
    # creating a reference of the default stdout
    old_stdout = sys.stdout

    # creating a container storing stdout
    gene_result = StringIO()
    sys.stdout = gene_result

    # calling ncbi gene database
    c = entrezpy.conduit.Conduit(email)
    fetch_genes = c.new_pipeline()
    sid = fetch_genes.add_search(
        {
            "db": "gene",
            "term": f"{accession}",
            "rettype": "count",
            "datetype": "pdat",
        }
    )
    fid = fetch_genes.add_fetch(
        {"retmax": 10, "retmode": "text", "rettype": "fasta"}, dependency=sid
    )
    c.run(fetch_genes)
    sleep(buffer)

    # store the string from the stdout into variable
    gene_conts = gene_result.getvalue()

    # now redirect back stdout to screen
    sys.stdout = old_stdout

    # parsing gene response contents
    gene_info = _parse_ncbi_genes_response(gene_conts)

    return gene_info


def _parse_ncbi_genes_response(contents: str) -> dict:
    """Parses ncbi's genes response

    Parameters
    ----------
    contents : str
        raw ncbi genes response

    Returns
    -------
    dict
        parsed gene meta data into a dictionary
    """
    split_contents = contents.splitlines()
    split_contents = [line_data for line_data in split_contents if len(line_data) != 0]
    chunked_contents = [
        split_contents[i : i + 6] for i in range(0, len(split_contents), 6)
    ]

    # TODO: need to redesign annotatation extraction is not general
    gene_dict = defaultdict(None)
    for idx, chunked_content in enumerate(chunked_contents):
        subdict = {}
        if "discontinued" in " ".join(chunked_content):
            continue

        if len(chunked_content) != 6:
            continue

        name = chunked_content[0].split(".")[-1].strip()
        gene_id = chunked_content[-1].split(":")[-1].strip()
        for line in chunked_content:

            # iterating through contents
            if line.startswith("Annotation"):
                if "complement" in line:
                    annotation = (
                        line.split(":")[-1]
                        .replace("(", "")
                        .replace(")", "")
                        .split()[1]
                        .replace(",", "")
                        .split("..")
                    )
                else:
                    annotation = (
                        line.split()[-1]
                        .replace("(", "")
                        .replace(")", "")
                        .replace(",", "")
                        .split("..")
                    )
                annotation_range = tuple([int(i) for i in annotation])
                subdict["annotation"] = annotation_range

            elif line.startswith("Other Designations:"):
                description = line.split(":")[-1].strip()
                subdict["name"] = f"{name}: {description}"

        gene_dict[gene_id] = subdict

    return gene_dict
