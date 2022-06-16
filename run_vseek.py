import json
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import dataframe_image as dfi

# VSeek imports
import vseek.common.io_files as vfiles
import vseek.common.loader as vloader
import vseek.common.vseek_paths as vsp
from vseek.utils.sequence_io import SequenceIO
from vseek.utils.sra_callers import download_fasta
from vseek.utils.vseek_analysis import dynamic_hamming
from vseek.utils.vseek_plots import plot_viral_composition, bat_country_geo_plot
from vseek.apis.ncbi import get_all_viral_accessions, get_viral_genes, get_viral_genomes


def rel_abundance(counts: int, counts_sum: int) -> float:
    """Calculates relative abundance given total sum

    Parameters
    ----------
    counts : int
        counts of virus
    counts_sum : int
        total amount of counts

    Returns
    -------
    float
        relative abundance percentage
    """
    return round((counts / counts_sum) * 100.0, 2)


if __name__ == "__main__":

    # CLI arguments
    parser = argparse.ArgumentParser(
        description="Command line program for characterizing bat viruses"
    )
    parser.add_argument(
        "-i", "--input", type=str, required=False, nargs="+", help="SRR id or profile if --profile flag is used"
    )
    parser.add_argument(
        "--profile", required=False, default=False, action="store_true", help="Path to viral counts json file"
    )
    parser.add_argument(
        "-e",
        "--email",
        default=None,
        type=str,
        required=False,
        help="Valid Email address to send requests to NCBI API. Note required if --test_run is flagged",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        default=0.40,
        type=float,
        required=False,
        help="Similarity thresh hold when identifying viruses",
    )
    parser.add_argument(
        "-rt",
        "--rel_threshold",
        default=1.0,
        type=float,
        required=False,
        help="relative abundance cutoff. The smaller the percentage, the noisier the data",
    )
    parser.add_argument(
        "--viral_counts",
        default=None,
        type=str,
        required=False,
        help="JSON file containing viral counts. Skips all downloading files",
    )
    args = parser.parse_args()

    # type checking
    if args.profile is False:
        if len(args.input) > 1:
            raise ValueError("Only 1 accession id can be provided")
        elif len(args.input) == 0:
            raise ValueError("Requires at least one accession number")

    elif args.profile is True:
        check = args.input[0]
        if not check.endswith(".json"):
            raise ValueError("Invalid file added. just be in json format")
        if len(args.input) > 1:
            raise ValueError("only one profile is allowed")
        args.input = args.input[0]

    if args.threshold > 1.0 or args.threshold <= 0:
        raise ValueError("Similarity threshold must be bewtween 0 <= x > 1.0")
    if args.rel_threshold > 100.0 or args.rel_threshold < 0:
        raise ValueError("Similarity threshold must be bewtween 0.0 < x > 100.0")


    # -----------------------
    # step 0. setup and data collection
    # -----------------------
    # initializing important folders
    database_path = vsp.init_db_path()
    results_path = vsp.init_results_dir()
    genome_dir = vsp.init_genome_db_path()
    fasta_dir = vsp.init_fasta_dir()
    prefetch_dir = vsp.init_prefetch_dir()
    profile_dir = vsp.init_profile_dir()
    string_db_dir = vsp.init_string_dir()

    # loading all datasets, if it doesn't exists, it will download it
    if args.email is not None and args.profile is False:
        print("Downloading all the datasets. Might take a while")
        download_fasta()
        get_viral_genomes()  # downloads all known viral genomes
        get_viral_genes(
            email=args.email, accession=args.accession
        )  # download viral genes with given genome

        viral_genes = vloader.load_viral_genes()
        meta_genome_paths = vfiles.get_meta_genomes_paths()
        ppi_df = vloader.load_human_ppi()  # loads human-viral database
        viral_accession_df = get_all_viral_accessions()  # downloading viral accessions
        dbat_vir_df = vloader.load_dbat_vir_db()
        ncbi_virus_data = vloader.load_bat_virus_data()
        geolocations_df = vloader.load_geolocations()
    else:
        # This will still download if files are actually missing
        dbat_vir_df = vloader.load_dbat_vir_db()
        ncbi_virus_data = vloader.load_bat_virus_data()
        geolocations_df = vloader.load_geolocations()
        genome_paths = vfiles.get_viral_genome_fasta_paths()
        meta_genome_paths = vfiles.get_meta_genomes_paths()
        ppi_df = vloader.load_human_ppi()  # loads human-viral database
        viral_accession_df = get_all_viral_accessions()  # downloading viral accessions

    # -----------------------
    # step 1. Discovery
    # -----------------------
    # if --test-run is true : skipping discovery step and use sample_data
    if args.profile is False:
        print(
            "WARNING: Running Discovery step, this will take a long time to finish. Rougly 6000y years for 39 million reads"
        )
        # parameters
        counts_save_path = (
            Path(results_path) / f"{args.accession}_viral_composition_counts.json"
        )
        all_accessions = list(vfiles.viral_genome_paths().keys())
        metagenome_path = vfiles.get_meta_genomes_paths()
        reader = SequenceIO(metagenome_path)
        reads = reader.lazy_load_fasta()
        threshold = (
            args.threshold
        )  # threshold between 40% and 80% is enough to get genius level

        counts = defaultdict(lambda: 0)
        for idx, read in enumerate(reads):

            if idx % 5 == 0:
                # saving every 5 analyzed reads
                with open(counts_save_path, "w") as outfile:
                    json.dump(counts, outfile)

            if idx % 100 == 0:
                print(f"Currently on read number {idx}")

            read_score = {}
            for acc_id in all_accessions:
                # load all the gene sequences
                gene_sequences = vloader.load_viral_genes(acc_id)

                top_score = 0.0
                for gene in gene_sequences:
                    score = dynamic_hamming(read=read.sequence, reference=gene)
                    if score >= threshold:
                        if score == 1:
                            top_score = score
                            break
                        elif score > top_score:
                            top_score = score

                read_score[acc_id] = top_score

            top_score_acc_ids = max(read_score, key=read_score.get)
            counts[top_score_acc_ids] += 1

        viral_count_data = vloader.load_viral_counts(counts_save_path)
    else:
        # NOTE: Place example counts here
        print("skipping Discovery step... loading")
        exmaple_counts_save_path = Path(database_path) / "viral_composition_counts.json"
        viral_counts_path = str(Path(args.input).absolute())
        viral_count_data = vloader.load_viral_counts(viral_counts_path)

    # -----------------------
    # step 2. Profiling
    # -----------------------
    geolocations_df = geolocations_df[["country", "name"]]
    geolocations_df.columns = ["iso_alpha", "country"]
    iso_codes_df = vloader.load_iban_iso_codes()

    # with open(viral_counts, "r") as infile:
    #     viral_count_data = json.load(infile)

    viral_count_df = pd.DataFrame.from_dict(
        data=viral_count_data, orient="index"
    ).reset_index()
    viral_count_df.columns = ["accession", "counts"]
    viral_count_sum = np.sum(viral_count_df["counts"])
    viral_count_df["rel_abundance"] = viral_count_df["counts"].apply(
        lambda count: rel_abundance(counts=count, counts_sum=viral_count_sum)
    )
    viral_count_df = viral_count_df.sort_values(by="rel_abundance", ascending=False)

    viral_count_profile_save = Path(profile_dir) / f"{args.input}"

    selected_counts = viral_count_df.loc[
        viral_count_df["rel_abundance"] >= args.rel_threshold
    ]
    viral_taxa = ncbi_virus_data[["Representative", "family", "genus"]]
    viral_taxa.columns = ["accession", "family", "genus"]
    identified_virus_df = viral_taxa.merge(selected_counts, on="accession")
    identified_virus_df = identified_virus_df.sort_values(
        by="rel_abundance", ascending=False
    )
    identified_virus_df = identified_virus_df.set_index("accession")

    # saving table as an image
    save_path_identified_virus_table_img = str(
        (Path(results_path) / "identified_virus_table.png").absolute()
    )

    save_path_identified_virus_table = str(
        (Path(results_path) / "identified_virus.csv").absolute()
    )
    identified_virus_df.to_csv(save_path_identified_virus_table, index=False)  # profile

    # PLOT 2: GEO Plot
    found_viral_fam = identified_virus_df["family"].unique().tolist()
    associated_bats_df = dbat_vir_df.loc[
        dbat_vir_df["Viral family"].isin(found_viral_fam)
    ]

    # removing any duplicates
    group = associated_bats_df.groupby("From Bat")

    list_dfs = []
    for name, df in group:
        new_df = df.drop_duplicates(subset=["Viral family", "From Bat"])
        list_dfs.append(new_df)

    associated_bats_df = pd.concat(list_dfs)
    associated_bats_df["Sampling country"].unique()

    # grouping by country and condensing data
    country_group = associated_bats_df.groupby("Sampling country")
    bat_country_data = []
    for name, df in country_group:

        if "congo" in name.lower():
            name = "Congo [Republic]"
        elif "united states of america" in name.lower():
            name = "United States"

        all_bats_str = " ".join(df["From Bat"].tolist())
        counts = len(df["From Bat"].tolist())
        viral_fam = ", ".join(df["Viral family"].unique().tolist())
        result = [name, viral_fam, all_bats_str, counts]

        bat_country_data.append(result)

    # constructing geo dataframe
    all_bats_country_df = pd.DataFrame(
        data=bat_country_data,
        columns=["country", "viral_fam", "all_bats", "bat_counts"],
    )
    all_bats_country_df = all_bats_country_df.sort_values(
        by="bat_counts", ascending=False
    )

    # augmenting all identified bats with location data
    alpha_iso_codes = iso_codes_df[["iso_alpha", "alpha_iso3"]]
    all_bats_country_df = all_bats_country_df.merge(geolocations_df, on="country")
    all_bats_country_df = all_bats_country_df.merge(alpha_iso_codes, on="iso_alpha")

    save_path_bats_country = Path(results_path) / "other_bats_country.csv"
    all_bats_country_df.to_csv(save_path_bats_country, index=False)  # geo profile

    # ASSOCIATED DISEASE PROFILE
    viral_disease_save_path = Path(results_path) / "viral_disease_profile.csv"
    ident_viral_acc = identified_virus_df.index.tolist()
    viral_info = ncbi_virus_data.loc[
        ncbi_virus_data["Representative"].isin(ident_viral_acc)
    ]
    viral_info.columns = ["accession"] + viral_info.columns.values.tolist()[1:]

    # merging relative counts and relative abundance
    viral_counts_df = identified_virus_df.reset_index()[["accession", "rel_abundance"]]
    viral_summary_df = viral_info.merge(viral_counts_df, on="accession").sort_values(
        "rel_abundance", ascending=False
    )
    viral_summary_df.columns = [
        "accession",
        "host",
        "disease_doc",
        "segment_name",
        "family",
        "genus",
        "taxon_id",
        "rel_abundance",
    ]
    viral_summary_df.to_csv(viral_disease_save_path)  # virus-disease profile

    # GENERATING AN INTERACTION PROFILE
    # filter the ppi data with the taxon id found in the viral_summary df
    found_virus_taxon_id = viral_summary_df["taxon_id"].unique().tolist()
    ppi_df = ppi_df.loc[ppi_df["species_2"].isin(found_virus_taxon_id)]
    vfiles.save_interaction_profiles(ppi_df)

    # -----------------------
    # step 3. Plotting
    # ----------------------
    viral_comp_save_path = (
        Path(results_path) / "viral_characterizations_piechart_plot.png"
    )
    plot_viral_composition(selected_counts, save_path=viral_comp_save_path)
    dfi.export(identified_virus_df, save_path_identified_virus_table_img)

    save_path_geo_plot = str((Path(results_path) / "geo_plot.png").absolute())
    bat_country_geo_plot(all_bats_country_df, save_path=save_path_geo_plot)

    print("process compelte")
