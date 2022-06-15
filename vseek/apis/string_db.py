import os
import gzip
from pathlib import Path

import requests
import pandas as pd
import numpy as np

import vseek.common.vseek_paths as vsp
import vseek.common.loader as vloader


def download_species_atlas() -> pd.DataFrame:
    """Downloads species atlas from string_db

    Returns
    -------
    pd.DataFrame
        species atlas
    """
    url = "http://viruses.string-db.org/download/species.v10.5.txt"
    resp = requests.get(url, timeout=10)
    resp.raise_for_status()

    species_df = pd.DataFrame([line.split("\t") for line in resp.text.splitlines()])
    header = ["taxon_id", "string_type", "string_name_compact", "official_name_ncbi"]
    species_df = species_df.iloc[1:]
    species_df.columns = header
    species_df

    # save the file in the database
    save_path = Path(vsp.init_string_dir()) / "stringdb_species_codes.tsv"
    species_df.to_csv(save_path, sep="\t", index=False)

    return species_df


def download_all_human_interactions() -> str:
    zipped_path = Path(vsp.ppi_db_path()) / "9606.protein.links.v10.5.txt.gz"
    ppi_save_path = Path(vsp.ppi_db_path()) / "human_ppi.tsv.gz"

    try:
        if not ppi_save_path.is_file():
            print("\nDownloading all human protein interactions with all species")
            url = "http://viruses.string-db.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz"
            resp = requests.get(url)
            resp.raise_for_status()

            with open(zipped_path, "wb") as outfile:
                outfile.write(resp.content)

            # removing content
            resp = ""
            del resp

            # unzipping file
            print("Constructing ppi DataFrame")
            conts = _ppi_contents_parser(zipped_path)

            # constructing ppi dataframe
            human_ppi_df = pd.DataFrame(
                data=conts,
                columns=["species_1", "protein_1", "species_2", "protein_2", "score"],
            )

            # reassign and deleting list
            row_data = ""
            del row_data

            # removing compressed file and raw_data
            os.remove(str(zipped_path.absolute()))

            # filtering down to only viral and protein interactions and annotation
            human_ppi_df = _filter_viral_interactions(human_ppi_df)

            human_ppi_df.to_csv(ppi_save_path, sep="\t", index=False)

            return human_ppi_df
    except EOFError:
        print("Incomplete download captured, downloading again ...")
        os.remove(str(zipped_path.absolute()))
        human_ppi_df = download_all_human_interactions()
        return human_ppi_df


def _filter_viral_interactions(ppi_interactions: pd.DataFrame) -> pd.DataFrame:
    """Filters to only protein and viral human interactions

    Returns
    -------
    _type_
        _description_
    """
    print("Filtering to only viral-human protein interactions")
    bat_vir_data = vloader.load_bat_virus_data()
    viral_taxon_id = bat_vir_data["taxon_id"].unique().astype(str)
    filtered_data = ppi_interactions.loc[
        ppi_interactions["species_2"].isin(viral_taxon_id)
    ]

    # annotate
    filtered_data = _annotate_viral_taxon(filtered_data)

    return filtered_data


def _annotate_viral_taxon(ppi_interactions: pd.DataFrame) -> pd.DataFrame:
    """Annotates viral information viral taxon id

    Parameters
    ----------
    ppi_interactions : pd.DataFrame
        viral-host interactions

    Returns
    -------
    pd.DataFrame
        annotated viral-host interactions
    """
    # annotation function used for the .apply() method (below)
    print("Annotating viral information")

    def _annotate_taxon(taxon_id: int, species_atlas: dict) -> str:
        """FOR PANDAS APPLY FUNCTION ONLY: annotation based on taxon id

        Parameters
        ----------
        taxon_id : int
            taxon id
        species_atlas : dict
            dictionary containing taxon id and annotation as key value pairs

        Returns
        -------
        str
            annotation
        """
        try:
            result = species_atlas[taxon_id]
            return result
        except KeyError:
            return np.nan

    # generating species_atlas for annotation
    species_anno_df = vloader.load_species_atlas()
    species_anno = dict(
        zip(species_anno_df["taxon_id"], species_anno_df["official_name_ncbi"])
    )

    # annotating with viral description
    # uses _annotate_taxon() function (above)
    ppi_interactions["annotation"] = (
        ppi_interactions["species_2"]
        .astype(int)
        .apply(lambda taxon_id: _annotate_taxon(taxon_id, species_atlas=species_anno))
    )
    return ppi_interactions


# loader module


def _ppi_contents_parser(ppi_file: str) -> list[str]:
    """Creates a generator of parsed ppi contents

    Parameters
    ----------
    ppi_file : str
        path pointing to compressed ppi file

    Returns
    -------
    list[str]
        contents of ppi file in a list
    """
    with gzip.open(ppi_file, "rb") as gz_file:

        for i, conts in enumerate(gz_file):

            # removing original headers
            if i == 0:
                continue

            # splitting contents
            row = conts.decode("utf-8").split()
            cont1 = row[0].split(".")
            cont2 = row[1].split(".")
            score = row[-1]
            new_row = cont1 + cont2 + [score]
            yield new_row
