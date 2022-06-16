# Read me

# VSeek: A human-bat virus characterization application.

`VSeek` is a simple command line applications that characterizes the viral composition from a meta-genomic sequence reads. It uses possess its own databases that uses data from [DBatVir](http://www.mgc.ac.cn/DBatVir/), [String.virus](http://viruses.string-db.org) and [NCBI genome](https://www.ncbi.nlm.nih.gov/genome/viruses/) data bases in order to generate in robust profiles that are useful for pandemic prevention. 

The program possess three mains steps:

- Disocvery → Deploys a searching algorithm to identify which viral
- Profiling → Uses the obtained viral counts from the Disocvery steps and generates profiles. (Viral-disease , geo  and protein-protein interaction profiles)
- Plotting → Uses the generated profiles to produce plots

## Dependencies

This program has only been tested on **Unix** based OS systems like **Apple** or **Linux**. Installations processes **may fail** it conducted in a Windows computer. 

If Windows is the only computer available to you, you can activate Microsoft windows system for Linux (WSL) in order to create a Linux based environment. Here is the [installation guide](https://docs.microsoft.com/en-us/windows/wsl/install). 

The dependencies and versions are listed in the `environment.yml` file. These dependencies are installed in this installation step. 

```yaml
name: "VSeek"
channels:
  - conda-forge
  - bioconda
  - anaconda
  - default
dependencies:
  - python>=3.10.0
  - bioconda::sra-tools==2.10.0
  - bioconda::entrezpy>=2.1.3
  - conda-forge::numpy>=1.22.0
  - conda-forge::matplotlib>=3.5.0
  - conda-forge::pandas>=1.4.0
  - conda-forge::seaborn>=0.11.0
  - conda-forge::lxml>=4.9.0
  - conda-forge::selenium>=4.2.0
  - conda-forge::webdriver-manager>=3.7.0
  - conda-forge::ipykernel>=6.13.0
  - conda-forge::dataframe_image>=0.1.0
  - conda-forge::plotly>=5.8.0
  - conda-forge::python-kaleido>=0.2.0
```

## Installation

This installation processes assumes that you already have the anaconda package management into your system system. If you do not, you can refer anaconda installation [page](https://docs.anaconda.com/anaconda/install/).  

First, clone this repository to your current directory by inputting: 

```
git clone https://github.com/axiomcura/VSeek.git
```

The next step is to create an environment and download all the dependencies that `VSeek` requires in order to conduct its analysis. `VSeek` has it’s own `environment.yml` file that contains the required programs and it’s version to download. this is executed by typing 

```
conda env create -f environment.yml
```

Then to activate the environment type:

```
conda activate VSeek
```

**Warning:** Your must be in this environment in order to run `VSeek`. The application will not work if the incorrect environment is loaded. 

To know if you are in the right environment, you should see `(VSeek)` your terminal behind the username. 

## Documentations

To access the documentation 

```
(VSeek) erikserrano@Eriks-MacBook-Pro-2 VSeek % python run_vseek.py --help
usage: run_vseek.py [-h] [-i INPUT [INPUT ...]] [--profile] [-e EMAIL] [-t THRESHOLD] [-rt REL_THRESHOLD] [--viral_counts VIRAL_COUNTS]

Command line program for characterizing bat viruses

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        SRR id or profile if --profile flag is used
  --profile             Path to viral counts json file
  -e EMAIL, --email EMAIL
                        Valid Email address to send requests to NCBI API. Note required if --test_run is flagged
  -t THRESHOLD, --threshold THRESHOLD
                        Similarity thresh hold when identifying viruses
  -rt REL_THRESHOLD, --rel_threshold REL_THRESHOLD
                        relative abundance cutoff. The smaller the percentage, the noisier the data
  --viral_counts VIRAL_COUNTS
                        JSON file containing viral counts. Skips all downloading files
```

## Use case

### Simple use case

This use case example avoids the usage of the Discovery steps as it takes an extensive amount of time to finish. In this use case example, we’ll be using the `--profile` parameter, which requires a viral-count profile in json format. 

```
python run_vseek.py -i SRR12464727_viral_composition_counts.json --profile
```

In your terminal, you should see these outputs:

```
skipping Discovery step... loading

Generating viral-human protein interaction plots
interactive plot saved: file:///Path/To/plots/Varicella-zoster_virus_%28strain_Dumas%29_ppi_plot.html
interactive plot saved: file:///Path/To/plots/Variola_virus_%28isolate_HumanIndiaInd31967%29_ppi_plot.html
interactive plot saved: file:///Path/To/plots//Human_herpesvirus_6A_%28strain_Uganda-1102%29_ppi_plot.html
interactive plot saved: file:///Path/To/plots/Human_cytomegalovirus_%28strain_Merlin%29_ppi_plot.html
Copy and paste these URI's into your browser

Creating protein-protein network profiles files in SIF and TXT formats

objc[55272]: Class WebSwapCGLLayer is implemented in both /System/Library/Frameworks/WebKit.framework/Versions/A/Frameworks/WebCore.framework/Versions/A/Frameworks/libANGLE-shared.dylib (0x7ffb5c7eb948) and /Applications/Google Chrome.app/Contents/Frameworks/Google Chrome Framework.framework/Versions/102.0.5005.115/Libraries/libGLESv2.dylib (0x1124cf290). One of the two will be used. Which one is undefined.
[0616/113616.889654:INFO:headless_shell.cc(660)] Written to file /var/folders/88/kxtrpbsj7mg68w7mlprrsbvr0000gn/T/tmphoymc4j9/temp.png.

Saving interactive geo plot
to view interactive plot: copy and pase this onto your browser:
file:///Path/To/plots/geo_plot.html
```

Two main stdouts to pay attention too is the `uri` link that is provided. This link allows you to see the interactive plot in your browser. This is done by copying the `uri` link and pasting it into your browser’s address bar.

### Produced plots and profiles

In the `./results` folder, will produces these files after running the command above

```
SRA_prefetch                              other_bats_country.csv
complete_ppi_profile.csv                  plots
fasta_files                               ppi-profiles
identified_virus.csv                      viral_disease_profile.csv
identified_virus_table.png                
```

- `SRA_prefetch` → stores the SRR prefetched files using `sra-toolk`'s `prefetch` executable
- `fasta_files` → stores sequences using `sra-toolkit`'s `fasterq-dump` executable
- `plots` → contains the are interactive and static images including interactive profiles, geo plots and viral composition plot
- `viral_characterizations_piechart_plot.png` is a pie chart displaying the diversity of identified viruses
- `identified_virus.csv` → are the identified virus profile and `identified_virus_table.png` is just a table image of the profile
- `viral_disease_profile.csv` → profile containing identified viruses and associated diseases
- `ppi-profiles` contains the interaction information in either txt or SIF (useful for [Cyotoscape](https://cytoscape.org)) to visualize protein-protein interactions

## Additional notes

- if you would like to conduct the discovery step, change the command line inputs to this:

```
python run_vseek.py SRR12464727
```

- you can build the VSeek’s database by using `build_database.py`

```
python run_vseek.py build_database.py
```

- If you are getting `ModuleNotFound` errors, conduct a full reinstall but removing he environment by typing `conda info --envs` and deleting the `VSeek`'s env path:
    
    ```
    (VSeek) User % conda info --envs
    # conda environments:
    #
    
    VSeek                 *  /Path/To/VSeek_env
    
    (VSeek) User % rm -rf /Path/To/VSeek_env
    ```
    
    - make sure not to put `Path/To/VSeek_env` in your terminal, this is used as a place holder for your real path when you execute `conda info --envs`
- Then you can do a complete reinstall by repeating the `installation` steps
