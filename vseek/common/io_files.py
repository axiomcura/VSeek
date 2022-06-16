import os
import glob
import json
from typing import Union
from collections import defaultdict
from pathlib import Path
import pandas as pd

# vseek imports
import vseek.common.vseek_paths as vsp
import vseek.common.loader as vloader
from vseek.common.errors import *
from vseek.common.checks import prefetch_dir_exists, metagenome_dir_exists


def genome_dir_paths() -> list[str]:
    """Obtains all viral genome directories in the genome database

    Returns
    -------
    list[str]
        list of paths to each fasta file in the genome database
    """
    genome_db_path = vsp.genome_db_path()
    query = f"{genome_db_path}/*"
    all_genome_dirs = glob.glob(query)

    return all_genome_dirs


def get_viral_genome_fasta_paths(query=None) -> dict:
    """Retrieves all genome fasta file paths

    query : str, list[str]
        takes in a single query (str) or multiple queries (list[str]).
        default query=None returns all paths

    Return
    ------
    dict
        accession id and viral genome fasta path as key value pairs
    """
    if query is None:
        return _viral_genome_fasta_path_lookup()
    else:
        fasta_file_paths = _viral_genome_fasta_path_lookup()
        if not isinstance(query, list):
            query = query.split()

        result_query = defaultdict(None)
        for q in query:
            result = fasta_file_paths[q]
            if result is None:
                raise FastaFileNotFound(f"{q} does not exist")
            result_query[q] = result

        return result_query


def get_genome_genes_paths(query=None) -> dict:
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
        return _genes_path_lookup()
    else:
        fasta_file_paths = _genes_path_lookup()
        if not isinstance(query, list):
            query = query.split()

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
        all_genome_profiles = _genome_dir_path_lookup()
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


def get_meta_genomes_paths() -> dict:
    """Gets all paths that points to downloaded meta-genome files

    Returns
    -------
    list[str]
        list of paths pointing to meta-genome files
    """
    check = metagenome_dir_exists()
    if check is False:
        raise FileNotFoundError("Meta-genome folder has not been created")

    # get all metagenome files
    all_metagenome_files = _all_metagenome_files()

    return all_metagenome_files


# -----------------------------
# Writer functions
# -----------------------------
def save_genome(accession: str, contents: str) -> None:
    """Saves genome into the genome database

    Parameters
    ----------
    accession : str
        accession number of the viral genome
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


def save_genes(accession: str, contents: dict) -> None:
    """Saves viral genome genes into a json file.

    Parameters
    ----------
    accession : str
        accession number of the viral genome
    contents : dict
        contains gene meta data

    Returns
    -------
    None
        saves a json file in the genome database
    """
    json_conts = defaultdict(None).fromkeys([accession])
    json_conts[accession] = contents

    save_path = Path(vsp.genome_db_path()) / accession / f"{accession}_genes.json"
    with open(save_path, "w") as infile:
        json.dump(json_conts, infile)


def save_interaction_profiles(ppi_df: pd.DataFrame):
    """Generates interaction profiles in SIF format

    Parameters
    ----------
    ppi_df : pd.DataFrame
        protein-protein interactions profile
    """
    print("\nCreating protein-protein network profiles files in SIF and TXT formats\n")
    s_atlas = vloader.load_species_atlas()
    x = ppi_df.groupby(by=["species_2"])
    for species_id, species_df in x:
        species_name = s_atlas.loc[s_atlas["taxon_id"] == species_id][
            "official_name_ncbi"
        ].iloc[0]
        species_name = "_".join(species_name.split()).replace("/", "")

        group2 = species_df.groupby("protein_2")

        interactions = []
        adjacency_list_conts = []
        for protein_name, interaction_df in group2:
            all_human_proteins = " ".join(interaction_df["protein_1"].tolist())
            sif_interaction_str = f"{str(protein_name)} pp {all_human_proteins}"
            interactions.append(sif_interaction_str)

        interaction_save_path = (
            Path(vsp.init_profile_dir()) / f"{species_name}_interaction.sif"
        )
        with open(interaction_save_path, "w") as outfile:
            for interaction in interactions:
                outfile.write(f"{interaction}\n")

        adjacency_save_path = (
            Path(vsp.init_profile_dir()) / f"{species_name}_interaction.txt"
        )
        with open(adjacency_save_path, "w") as outfile:
            for adj_interaction in adjacency_list_conts:
                outfile.write(f"{adj_interaction}\n")


# -----------------------------
# Removers
# -----------------------------
def clean_all_genes() -> None:
    """Removes all viral gene profile files

    Returns
    --------
    None
    """
    targets = get_genome_genes_paths().values()
    if len(targets) == 0:
        print("Warning: There are no gene files")
        return

    for target in targets:
        if Path(target).is_file():
            os.remove(target)
        else:
            print(f"Warning: {target} is not a file")


# -----------------------------
# Private functions
# -----------------------------
def _viral_genome_fasta_path_lookup() -> dict:
    """Creates a dictionary of all fasta paths associated with accession number

    Returns
    -------
    dict
        accession number and path as key value pairs
    """
    genome_db_paths = genome_dir_paths()

    all_fasta_paths = defaultdict(None)
    for genome_path in genome_db_paths:
        genome_id = genome_path.split("/")[-1].split(".")[0]
        query = f"{genome_path}/*.fasta"
        fasta_file_path = glob.glob(query)[0]
        all_fasta_paths[genome_id] = fasta_file_path

    return all_fasta_paths


def _genes_path_lookup() -> dict:
    """Creates a dictionary of all genome profile paths associated with accession number

    Returns
    -------
    dict
        accession number and profile path as key value pairs
    """
    genome_db_paths = genome_dir_paths()

    all_genome_profile_paths = defaultdict(None)
    for genome_path in genome_db_paths:
        genome_id = genome_path.split("/")[-1]
        query = f"{genome_path}/*_genes.json"
        fasta_file_path = glob.glob(query)[0]

        if len(fasta_file_path) == 0:
            continue
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


def _all_metagenome_files():
    metagenome_dir = vsp.metagenome_path()
    query = f"{metagenome_dir}/*.fasta"
    all_files = glob.glob(query)
    if len(all_files) == 1:
        return all_files[0]

    return all_files
