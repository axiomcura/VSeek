from pathlib import Path
import requests
import pandas as pd

# VSeek imports
import vseek.common.vseek_paths as vsp
from vseek.utils.parsers import parse_ncbi_viral_accessions


def get_all_viral_accessions() -> pd.DataFrame:
    """Sends a request to ncbi to obtain all viral genome accessions
    """

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




