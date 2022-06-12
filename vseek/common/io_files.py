import glob
from collections import defaultdict
from pathlib import Path

# vseek imports
import vseek.common.vseek_paths as vsp
from vseek.common.errors import *
from vseek.common.checks import prefetch_dir_exists


def genome_dir_paths() -> list[str]:
    """Obtains all genome directories in the genome database

    Returns
    -------
    list[str]
        list of paths to each fasta file in the genome database
    """
    genome_db_path = vsp.genome_db_path()
    query = f"{genome_db_path}/*"
    all_genome_dirs = glob.glob(query)

    return all_genome_dirs


def get_genome_fasta_paths(query=None) -> dict:
    """Retrieves all genome fasta file paths

    query : str, list[str]
        takes in a single query (str) or multiple queries (list[str]).
        default query=None returns all paths

    Return
    ------
    dict
        accession id and genome profile paths as key value pairs
    """
    if query is None:
        return _fasta_path_lookup()
    else:
        fasta_file_paths = _fasta_path_lookup()
        if not isinstance(query, list):
            query = list(query)

        result_query = defaultdict(None)
        for q in query:
            result = fasta_file_paths[q]
            if result is None:
                raise FastaFileNotFound(f"{q} does not exist")
            result_query[q] = result

        return result_query


def get_genome_profile_paths(query=None) -> dict:
    """Retrieves all genome profile file paths

    query : str, list[str]
        takes in a single query (str) or multiple queries (list[str]).
        default query=None returns all paths

    Return
    ------
    dict
        accession id and genome profile paths as key value pairs
    """
    if query is None:
        return _profile_path_lookup()
    else:
        fasta_file_paths = _fasta_path_lookup()
        if not isinstance(query, list):
            query = list(query)

        result_query = defaultdict(None)
        for q in query:
            result = fasta_file_paths[q]
            if result is None:
                raise ProfileNotFound(f"{q} does not exist")
            result_query[q] = result

        return result_query


def get_genome_dir_path(query=None) -> dict:
    """Gets the genome profiles
    query : str, list[str]
        takes in a single query (str) or multiple queries (list[str]).
        default query=None returns all paths

    Returns
    -------
    list[str], str
        returns path with given query or returns a list of paths
        if query=None, all paths will be returned
    """
    if query is None:
        return _genome_dir_path_lookup()
    else:
        all_genome_profiles = _fasta_path_lookup()
        if not isinstance(query, list):
            query = list(query)

        result_query = defaultdict(None)
        for q in query:
            result = all_genome_profiles[q]
            if result is None:
                raise GenomeDirectoryNotFound(f"{q} does not exist")
            result_query[q] = result

        return result_query


def get_prefetch_files() -> list[str]:
    """Gets all paths that points to downloaded prefetch files

    Returns
    -------
    list[str]
        list of paths pointing to sra prefetched files
    """
    # checking the prefetch folder exists
    check = prefetch_dir_exists()
    if check is False:
        raise FileNotFoundError("Prefetch folder has not been created")

    # get all prefetch files
    all_pfiles = _all_prefetch_files()

    return all_pfiles


# -----------------------------
# Writer functions
# -----------------------------
def save_genome(accession: str, contents: str) -> None:
    """Saves genome into the genome database

    Parameters
    ----------
    accession : str
        accession number of the genome
    contents : str
        Viral genome Fasta contents

    Return
    ------
    None
        Saves files into genome database "./db"
    """
    # writing out fasta file
    genome_path = Path(vsp.init_genome_db_path()) / accession
    genome_path.mkdir(exist_ok=True)

    save_path = genome_path / f"{accession}.fasta"
    with open(save_path, "w") as outfile:
        outfile.write(contents)


# -----------------------------
# Private functions
# -----------------------------
def _fasta_path_lookup() -> dict:
    """Creates a dictionary of all fasta paths associated with accession number

    Returns
    -------
    dict
        accession number and path as key value pairs
    """
    genome_db_paths = genome_dir_paths()

    all_fasta_paths = defaultdict(None)
    for genome_path in genome_db_paths:
        genome_id = genome_path.split("/")[-2]
        query = f"{genome_path}/*.fasta"
        fasta_file_path = glob.glob(query)
        all_fasta_paths[genome_id] = fasta_file_path

    return all_fasta_paths


def _profile_path_lookup() -> dict:
    """Creates a dictionary of all genome profile paths associated with accession number

    Returns
    -------
    dict
        accession number and profile path as key value pairs
    """
    genome_db_paths = genome_dir_paths()

    all_genome_profile_paths = defaultdict(None)
    for genome_path in genome_db_paths:
        genome_id = genome_path.split("/")[-2]
        query = f"{genome_path}/*.json"
        fasta_file_path = glob.glob(query)
        all_genome_profile_paths[genome_id] = fasta_file_path

    return all_genome_profile_paths


def _genome_dir_path_lookup() -> dict:
    """Creates a dictionary of all genome directory paths associated with accession number

    Returns
    -------
    dict
        accession number and genome path as key value pairs
    """
    all_profiles = defaultdict(None)
    for gdb_file_path in genome_dir_paths():
        genome_id = gdb_file_path.rsplit("/", 1)[-1]
        all_profiles[genome_id] = gdb_file_path

    return all_profiles


def _all_prefetch_files() -> list[str]:
    """Returns a list of files sra files inside the prefetch directory

    Returns
    -------
    list[str]
        list of paths pointing to prefetch files
    """
    prefetch_dir = vsp.prefetch_path()
    query = f"{prefetch_dir}/*/*.sra"
    all_files = glob.glob(query)

    return all_files
