import os
from pathlib import Path


def root_path() -> str:
    """Returns project root path

    Returns
    -------
    str
        project root path
    """
    splitter = "VSeek"
    working_dir = os.getcwd()
    root_path = working_dir.split(splitter)[0] + splitter

    return root_path


def db_path() -> str:
    """ Returns path pointing to database folder
    
    Returns
    -------
    str
        path to db directory
    """
    root = Path(root_path())
    db_path = str(root / "db")

    return db_path
    

def init_db_path() -> str:
    db_path_obj = Path(db_path())
    db_path_obj.mkdir(exist_ok=True)

    return str(db_path_obj.absolute())
