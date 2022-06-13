def parse_ncbi_viral_accessions(contents: str) -> tuple:
    """Parses all viral accession

    Parameters
    ----------
    contents : str
        single string obtained from the request object. Contains all the information
        in

    Returns
    -------
    tuple
        header and clean list of viral accession information.
    """
    contents = contents.splitlines()
    cols = contents[1].strip().replace("\t", "-").replace('"', '').split("-")[1:]
    all_data = []
    for line_cont in contents[2:]:
        data = line_cont.strip().replace("\n", "").replace("\t", "---").replace('"', '').split("---")
        all_data.append(data)


    return (cols, all_data)


def parse_ncbi_genes_response(contents: str) -> None:
    """Parses ncbi's genes response

    Parameters
    ----------
    contents : str
        raw ncbi genes response
    """