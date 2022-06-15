from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px

#VSeek imports
import vseek.common.vseek_paths as vsp

def plot_viral_composition(viral_counts_df: pd.DataFrame, save_path: str) -> None:
    """Creates a viral composition pie chart

    Parameters
    ----------
    viral_counts_df : pd.DataFrame
        contains information of viral composition from sample
    save_path : str
        path to save image
    """
    labels = viral_counts_df["accession"]
    scores = viral_counts_df["counts"]

    # generate a pie chart
    colors = sns.color_palette('pastel')[0:10]

    #create pie chart
    plot = plt.pie(scores, labels=labels, colors=colors, autopct='%.0f%%', textprops={"fontsize": 8})
    plt.title("Viral Characterizations")
    plt.savefig(save_path)


def bat_country_geo_plot(bat_country_df: pd.DataFrame, save_path: str) -> None:
    """Generate a geo plot where all bats are located with identified virus.
    Also generates an interactive plot with .html

    Parameters
    ----------
    bat_country_df : pd.DataFrame
        data frame that contains identified virus and their locations
    save_path : str
        path to save image

    """
    fig = px.scatter_geo(bat_country_df, locations="alpha_iso3" , color="bat_counts",
                        title="Identified Bat's Viral Family strains and Location",
                        hover_name="country", size="bat_counts",
                        projection="natural earth", hover_data=[bat_country_df["country"], bat_country_df["viral_fam"]])

    fig.update_layout(height=1000, margin={"r":1,"t":100,"l":30,"b":0}, font=dict(size=25))

    # saving as interactive plot
    save_interactive_path = save_path.replace(".png", ".html")
    fig.write_html(save_interactive_path)

    # saving original plot
    with open(save_path, "wb") as outfile:
        img_bytes = fig.to_image(format="png", engine="kaleido")
        outfile.write(img_bytes)
