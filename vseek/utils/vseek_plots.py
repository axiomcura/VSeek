from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import networkx
from bokeh.io import output_file, show, save
from bokeh.models import Range1d, Circle, ColumnDataSource, MultiLine
from bokeh.plotting import figure
from bokeh.plotting import from_networkx

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
    print("\nSaving interactive geo plot")
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

    print(f"to view interactive plot: copy and pase this onto your browser:\n{Path(save_interactive_path).as_uri()}")

def ppi_interactive_plot(name: str, ppi_df: pd.DataFrame) -> None:
  """Generates an interactive plot that can be viewed in the browser.

  Parameters
  ----------
  name : str
      Name of the species
  ppi_df : pd.DataFrame
      pd.DataFrame that contains edge connections [viral_protein, human_protein]
      as entries
  """
  title = f'{name} and human protein-protein interactions network'


  # create edge graph
  G = networkx.from_pandas_edgelist(ppi_df, "protein_2", "protein_1")

  #Create a plot — set dimensions, toolbar, and title
  HOVER_TOOLTIPS = [("Protein_id", "@index")]
  plot = figure(tooltips = HOVER_TOOLTIPS,
                tools="pan,wheel_zoom,save,reset", active_scroll='wheel_zoom',
              x_range=Range1d(-10.1, 10.1), y_range=Range1d(-10.1, 10.1), title=title,
              plot_width=1300, plot_height=800)

  network_graph = from_networkx(G, networkx.spring_layout, scale=10, center=(0, 0))

  #Set node size and color
  network_graph.node_renderer.glyph = Circle(size=15, fill_color='skyblue')

  #Set edge opacity and width
  network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)

  #Add network graph to the plot
  plot.renderers.append(network_graph)


  save_path = Path(vsp.init_plots_dir()) / f"{name}_ppi_plot.html"
  output_file(str(save_path.absolute()))
  print(f"interactive plot saved: {save_path.as_uri()}")
  save(plot)