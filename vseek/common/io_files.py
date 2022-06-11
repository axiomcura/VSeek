import glob
from pathlib import Path

# vseek imports
import vseek.common.vseek_paths as vsp

def genome_db_files() -> list[str]:
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