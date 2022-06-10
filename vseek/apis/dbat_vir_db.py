import os
from pathlib import Path
import shutil
import tempfile
import glob
import warnings
from time import sleep
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

# VSeek imports
import vseek.common.vseek_paths as vsp


# ignoring deprecation warnings
warnings.filterwarnings("ignore")


def _group_dfs(data_path: str) -> list:
    """Places all pandas data frames into a list

    Parameters
    ----------
    data_path : str
        path containing all DBatVir data

    Returns
    -------
    list
        list of pd.DataFrames of DBatVir data
    """
    print("Merging downloaded data into one .csv file")

    # check if all the files are downloaded (there should be a total of 6)
    # -- if not all were downloaded, it will print a warning
    all_files = glob.glob(f"{data_path}/*")
    if not len(all_files) == 6:
        print("WARNING: not all files were downloaded")

    # now for the pandas part (merging all downloaded data into a giant .csv file)
    # -- creating all
    df_list = []
    for data_file in all_files:
        df = pd.read_excel(data_file)

        # editing header and reformating dataset
        header = df.iloc[0]
        df = df[1:]
        df.columns = header
        df_list.append(df)

    return df_list


def collector() -> list:
    """Downloads all Bat associated viruses from DBatVir

    Returns
    -------
    list
       list of pd.DataFrames containing viral and bat data
    """

    # special options for chrome drivers
    # -- default options
    options = webdriver.ChromeOptions()
    options.add_argument("--window-size=1920,1080")
    options.add_argument("--disable-extensions")
    options.add_argument("--proxy-server='direct://'")
    options.add_argument("--proxy-bypass-list=*")
    options.add_argument("--start-maximized")
    options.add_argument("--headless")
    options.add_argument("--disable-gpu")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--no-sandbox")
    options.add_argument("--ignore-ssl-error=yes")
    options.add_argument("--ignore-certificate-errors")

    # -- we should also specify that we want to download the files to the current directory
    temp_dl_path = tempfile.TemporaryDirectory()
    dl_path = temp_dl_path.name

    prefs = {"download.default_directory": dl_path}
    options.add_experimental_option("prefs", prefs)

    # now let connect to the database (deploying bot)
    # -- with our options
    print("Entering Bat-associated Virus database")
    driver = webdriver.Chrome(
        service=Service(ChromeDriverManager().install()), chrome_options=options
    )
    driver.get("http://www.mgc.ac.cn/cgi-bin/ZOVER/mainTable.cgi?db=chiroptera")
    sleep(3)  # buffer time to load the page
    print("Connected to DBatVir!")

    # ------------------------------
    # In the website now
    # ------------------------------
    # It turns out that these clickable elements have unique IDs! (Easy to find)
    # -- we can tell the web bot to click on those
    # -- we are downloading per virus group and not per virus

    # group name and element id as key value pairs
    virus_groups = {
        "Retro-transcribing-viruses": "35268_anchor",
        "dsDNA-viruses-noRNA-stage": "35237_anchor",
        "dsRNA-virus": "35325_anchor",
        "ssDNA-virus": "29258_anchor",
        "ssRNA-negative-strand-viruses": "35301_anchor",
        "ssRNA-positive-strand-viruses-noDNA-Stage": "35278_anchor",
    }

    # now lets iterate!
    print("Downloading Bat-Virus data files")
    for virus_group, elm_id in virus_groups.items():

        # find group button and click
        group_btn = driver.find_element_by_id(elm_id)
        group_btn.click()
        sleep(3)

        # serching for download button and click
        print(f"Downloading: {virus_group} data...")
        dl_btn = driver.find_element_by_xpath(
            "//*[contains(text(), 'Save this table')]"
        )

        sleep(2)  # for every grouped action, add a sleep to prevent server overload
        dl_btn.click()
        sleep(1)

    # ending the bot life
    print("Download complete! Closing portal...")
    sleep(4)
    driver.quit()

    # getting list of and removing temporary directory
    dfs = _group_dfs(dl_path)
    shutil.rmtree(dl_path)

    return dfs


def collect_dbatvir_data() -> pd.DataFrame:
    """Main wrapper to download data from DBatVir"""

    # init db path
    db_path_str = vsp.db_path()
    db_path = Path(db_path_str)
    db_path.mkdir(exist_ok=True)
    save_path = os.path.join(db_path, "DBatVir_db.csv.gz")

    # if it exists in the database, no need to download
    check = Path(save_path).is_file()
    if check is True:
        return pd.read_csv(save_path)

    else:
        # Collect all data
        df_list = collector()

        # save data to db
        main_df = pd.concat(df_list, axis=0)
        main_df.to_csv(save_path, compression="gzip", index=False)

        return main_df
