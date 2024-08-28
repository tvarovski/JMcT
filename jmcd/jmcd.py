# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import logging
from typing import Literal
from numpy import NaN
from natsort import index_natsorted
from jmcd.settings import settings
from jmcd.settings import chromosomes
from jmcd.settings import am1003_to_roman, am1003_to_numeric, roman_to_numeric, numeric_to_roman, S288C_to_roman, S288C_to_numeric
from jmcd.settings import chr_centromeres

def getGenomeSize() -> int:
    """
    Calculates the size of the genome by adding up the lengths of all chromosomes.

    Returns:
        int: The total size of the genome.
    """

    #calculate genome size by adding up the lengths of all chromosomes
    genome_size = 0
    for chromosome in chromosomes:
        genome_size += chromosome[1]
    
    return genome_size

def load_and_prep_data(path: str) -> pd.DataFrame:
    """
    This function loads a CSV file from the given path, filters the data, and keeps only the necessary columns.

    Parameters:
    path (str): The path to the CSV file.

    Returns:
    pd.DataFrame: The processed DataFrame.

    """
    df =  pd.read_csv(path)
    df = df[df["Type"] == "SNV"]
    df = df[df["Reference allele"] == "No"]
    df = df[df["Chromosome"] != "MT_AM_1003_v1"]

    #keep only columns that are needed
    columns_keep =[
        'Chromosome',
        'Region',
        'Type',
        'Reference',
        'Allele',
        'Zygosity',
        'Count',
        'Coverage',
        'Frequency',
        'Probability',
        'Forward read count',
        'Reverse read count',
        'Average quality',
        'BaseQRankSum'
    ]
    df = df[columns_keep]

    logging.info(f"Loaded {path[:12]}... with {len(df)} unique SNVs.")

    return df

def csvDataFrameImport(data_path: str, sample_name_len: int = 8) -> list:
    """
    This function imports CSV files from a given directory, processes them, and appends them to a list. 
    Each DataFrame in the list represents a sample, with the sample name derived from the filename.

    Parameters:
    data_path (str): The path to the directory containing the CSV files.
    sample_name_len (int): The length of the sample name to be extracted from the filename. Default is 8. If set to 0, the entire filename will be used as the sample name.

    Returns:
    list: A list of DataFrames, each representing a sample.

    """
    frames_list_samples = []

    for filename in os.listdir(data_path):

        #if file is a directory, skip it
        if os.path.isdir(f"{data_path}/{filename}"):
            logging.info(f"'{filename}' is a directory. Skipping.")
            continue

        assert filename.endswith(".csv"), f"{filename} is not a csv file. Please make sure that all files in the {data_path} directory are csv files"

        if sample_name_len == 0:
            sample_name = filename
        else:
            sample_name = filename[:sample_name_len]
        
        logging.debug(f"Loading {filename} as sample {sample_name}.")
        sample_df = load_and_prep_data(f"{data_path}/{filename}")
        sample_df["Sample"] = sample_name

        frames_list_samples.append((sample_df, sample_name))

    logging.info(f"found {len(frames_list_samples)} samples in the {data_path} directory")
    logging.info(f"sample names: {[sample[1] for sample in frames_list_samples]}")

    return(frames_list_samples)

def findSpectraStrandwise(colA: str, colB: str) -> str:
    """
    Given two nucleotide bases, returns the type of substitution that occurs when
    transitioning from the first base to the second base.

    Args:
        colA (str): The first nucleotide base.
        colB (str): The second nucleotide base.

    Returns:
        str: The type of substitution that occurs when transitioning from colA to colB.
    """
    if ((colA == 'C') & (colB == 'T')):
        return('C_to_T')
    elif ((colA == 'C') & (colB == 'A')):
        return('C_to_A')
    elif ((colA == 'C') & (colB == 'G')):
        return('C_to_G')
    elif ((colA == 'T') & (colB == 'C')):
        return('T_to_C')
    elif ((colA == 'T') & (colB == 'G')):
        return('T_to_G')
    elif ((colA == 'T') & (colB == 'A')):
        return('T_to_A')
    
    elif ((colA == 'G') & (colB == 'A')):
        return('G_to_A')
    elif ((colA == 'G') & (colB == 'T')):
        return('G_to_T')
    elif ((colA == 'G') & (colB == 'C')):
        return('G_to_C')
    elif ((colA == 'A') & (colB == 'G')):
        return('A_to_G')
    elif ((colA == 'A') & (colB == 'C')):
        return('A_to_C')
    elif ((colA == 'A') & (colB == 'T')):
        return('A_to_T')

def renameChr(
        df: pd.DataFrame, 
        col_name: str, 
        rename_scheme: Literal[
            "am_to_roman", 
            "am_to_numeric",
            "roman_to_numeric",
            "numeric_to_roman",
            "S288C_to_roman",
            "S288C_to_numeric"],
        custom_scheme: dict[str,str] | None = None
        ) -> pd.DataFrame:
    """
    Given a DataFrame and a column name, renames the chromosome names in the column to either roman numerals or S288C format.

    Args:

        df (pd.DataFrame): The DataFrame containing the column to be renamed.
        col_name (str): The name of the column to be renamed.
        rename_scheme (str): The format to which the chromosome names should be renamed.
            options:
            "am_to_roman" - renames the chromosome names from AM1003 to roman numerals
            "am_to_numeric" - renames the chromosome names from AM1003 to S288C format
            "roman_to_numeric" - renames the chromosome names from roman numerals to S288C format
            "numeric_to_roman" - renames the chromosome names from S288C format to roman numerals
        custom_scheme (dict): A custom dictionary that can be used to rename the chromosome names formatted as {"old_name": "new_name"}.

    Returns:
        pd.DataFrame: The DataFrame with the renamed column.
    """

    if rename_scheme == "am_to_roman":
        rename_dict = am1003_to_roman
    elif rename_scheme == "am_to_numeric":
        rename_dict = am1003_to_numeric
    elif rename_scheme == "roman_to_numeric":
        rename_dict = roman_to_numeric
    elif rename_scheme == "numeric_to_roman":
        rename_dict = numeric_to_roman
    elif rename_scheme == "S288C_to_roman":
        rename_dict = S288C_to_roman
    elif rename_scheme == "S288C_to_numeric":
        rename_dict = S288C_to_numeric

    if custom_scheme is not None:
        #overwrite the default rename_dict with the custom_scheme
        logging.info(f"Using custom rename scheme: {custom_scheme}")
        rename_dict = custom_scheme
                
    df[col_name] = df[col_name].replace(rename_dict)
    return df

def add_clusters_patch_to_ax(ax, pval_list, pval_dict, sample_name):

    #assign a descrete colormap for pvalues from the seaborn's colormap
    #generate a list of colors from colormap
    #cmap = sns.color_palette("viridis", as_cmap=True)
    cmap = sns.color_palette("ch:s=-.2,r=.6", as_cmap=True).reversed()
    cmap_colors = cmap(np.linspace(0, 1, len(pval_list)))
    cmap_dict = {pval_list[i]: cmap_colors[i] for i in range(len(pval_list))}

    pval_counter = 0
    for pval in pval_list:

        cluster_df = pval_dict[pval].copy()
        
        cluster_df["Sample_Cluster_ID"] = cluster_df.index

        #convert the chromosome names in the cluster_df["Chromosome"] column to numeric
        cluster_df = renameChr(cluster_df, "Chromosome", "am_to_numeric")

        cluster_dict = {}
        for row in cluster_df.itertuples(index=False):
            Cluster_ID = row.Sample_Cluster_ID
            Chromosome = row.Chromosome
            Start_position = row.Start
            Cluster_Length = row.Length
            
            cluster_dict[Cluster_ID] = [Chromosome, Start_position, Cluster_Length]

        for key in cluster_dict:
            cluster_chromosome = cluster_dict[key][0]

            if not isinstance(cluster_chromosome, str):
                logging.warning(f"Cluster chromosome is not a string. Skipping cluster {key}. for sample {sample_name}.")
                continue
            
            cluster_chromosome = cluster_chromosome.replace("chr", "")
            cluster_chromosome = int(cluster_chromosome)
            cluster_chromosome_coord = cluster_chromosome - (1.45) + (pval_counter*0.04) #(1.5 - offset to move rectangle lower/higher)
            cluster_start = cluster_dict[key][1]
            cluster_length = cluster_dict[key][2]
            
            logging.debug(f"Drawing cluster {key} with pvalue {pval} and color {cmap_dict[pval]} at {cluster_chromosome}:{cluster_start}-{cluster_start+cluster_length}.")
            ax.add_patch(plt.Rectangle((cluster_start+500, cluster_chromosome_coord), cluster_length+500, 0.5, edgecolor='black', linewidth=0.1, fill=True, facecolor=cmap_dict[pval], alpha=1, zorder=-2))
        pval_counter += 1

def drawSNPMap(
        df_SNVs: pd.DataFrame,
        pval_dict: dict,
        sample_name: str,
        saveMap: bool = True,
        showMap: bool = False,
        ) -> None:
    
    #To rename legend elements, change plot markers/colors, modify here
    label_dict = {
        "cen": 'Centromere',
        "SPECTRA_STRANDWISE": "Reference Allele",
        "Position": "Position",
        "G": "G->N",
        "C": "C->N",
        "A": "A->N",
        "T": "T->N",
        "C_to_T": "C->T",
        "C_to_A": "C/G->",
        "C_to_G": "C/G->",

        "G_to_A": "G->A",
        "G_to_C": "C/G->",
        "G_to_T": "C/G->",

        'T_to_A': "A/T->N",
        'T_to_C': "A/T->N",
        'T_to_G': "A/T->N",
        'A_to_C': "A/T->N",
        'A_to_G': "A/T->N",
        'A_to_T': "A/T->N",
        }  
    markers = {
        'cen': "o",
        "C_to_T": "$|$",
        "C_to_A": "$|$",
        "C_to_G": "$|$",

        "G_to_A": "$|$",
        "G_to_C": "$|$",
        "G_to_T": "$|$",

        'T_to_A': "$|$",
        'T_to_C': "$|$",
        'T_to_G': "$|$",
        'A_to_C': "$|$",
        'A_to_G': "$|$",
        'A_to_T': "$|$",
        }
    palette = {
        'cen': "white",
        "C_to_T": "tab:red",
        "C_to_A": "tab:olive",
        "C_to_G": "tab:olive",

        "G_to_A": "tab:blue",
        "G_to_C": "tab:olive",
        "G_to_T": "tab:olive",

        'A_to_C': "tab:gray",
        'A_to_G': "tab:gray",
        'A_to_T': "tab:gray",

        'T_to_A': "tab:gray",
        'T_to_C': "tab:gray",
        'T_to_G': "tab:gray"
        }
    legend_order = [
        "cen", 
        "Position", 
        "C_to_T", 
        "G_to_A", 
        "C_to_A", 
        "C_to_G", 
        "G_to_C", 
        "G_to_T", 
        "T_to_A", 
        "T_to_C", 
        "T_to_G", 
        "A_to_C", 
        "A_to_G", 
        "A_to_T"]

    try:
        df_SNVs["SPECTRA_STRANDWISE"] = df_SNVs.apply(lambda x: findSpectraStrandwise(x["Reference"], x["Allele"]), axis=1)
        #calculate fraction a3a
        total_mutations = len(df_SNVs)
        a3a_mutations = len(df_SNVs[df_SNVs["SPECTRA_STRANDWISE"].isin(["C_to_T", "G_to_A"])])
        a3a_fraction = a3a_mutations / total_mutations
        logging.info(f"Sample: {sample_name}. Found {a3a_mutations} A3A spectra mutations out of {total_mutations} total SNVs. A3A fraction: {a3a_fraction}.")
        logging.debug(f'\n{df_SNVs["SPECTRA_STRANDWISE"].value_counts()}')
    except:
        logging.error(f"Error in findSpectraStrandwise function. Skipping sample {sample_name}. The sample dataframe is might be empty (size: {len(df_SNVs)}).")
        plt.close()
        return
    
    sns.set_style("ticks")
    sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})
    fig, ax = plt.subplots()

    from jmcd.settings import chromosomes
    df_chr_lengths = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])
    df_chr_lengths = renameChr(df_chr_lengths, "chromosome", "S288C_to_roman")
    df_chr_lengths.plot(kind='barh',legend=False, ax=ax, x="chromosome", y="end_position", fontsize=15, figsize=(20,6), edgecolor="black", linewidth=1, color="#E6E6E6" , zorder=2, title=sample_name, label="Position")

    df_SNVs = renameChr(df_SNVs, "Chromosome", "am_to_numeric")
    df_SNVs = pd.concat([df_SNVs, chr_centromeres], ignore_index=True)
    df_SNVs = df_SNVs.sort_values(by="Chromosome", key=lambda x: np.argsort(index_natsorted(df_SNVs["Chromosome"])))
    df_SNVs = renameChr(df_SNVs, "Chromosome", "numeric_to_roman")

    sns.scatterplot(ax=ax, data=df_SNVs[df_SNVs["Reference"].isin(["cen"])], x="Region", y="Chromosome", hue="SPECTRA_STRANDWISE", palette=palette, style="SPECTRA_STRANDWISE", markers=markers, s=100, alpha=1,zorder=2, linewidth=0.75, edgecolor="black")
    sns.scatterplot(ax=ax, data=df_SNVs[~df_SNVs["Reference"].isin(["cen"])], x="Region", y="Chromosome", hue="SPECTRA_STRANDWISE", palette=palette, style="SPECTRA_STRANDWISE", markers=markers, s=100, alpha=0.85,zorder=3, linewidth=0.05, edgecolor="black")

    pval_list = settings['pval_list']
    #create a hard copy of the pval_dict to avoid modifying the original dictionary and reverse the pval_list
    pval_list_plot = pval_list.copy()
    pval_list_plot.reverse()
    add_clusters_patch_to_ax(ax, pval_list_plot, pval_dict, sample_name)

    ax.set_xlabel("POSITION")
    ax.set_ylabel("CHROMOSOME")

    #dont repeat the same item in the legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))

    #order the legend
    while True:
        try:
            by_label = {label_dict.get(key): by_label[key] for key in legend_order}
            break

        except Exception as e:
            try:
                logging.warning(f"Could not relabel legend. {e}. Removing {e} from legend_order and trying again.")
                #extract the string from the error message
                e = str(e).split("'")[1]
                legend_order = [x for x in legend_order if x != str(e)]
            except Exception as e:
                logging.warning(f"Cound not relabel legend. {e}.")
                break

    ax.legend(by_label.values(), by_label.keys(), fontsize=10, loc='lower right', bbox_to_anchor=(0.95, 0.1), title=None)
    plt.tight_layout()

    if saveMap:
        plt.savefig(f'figures/{sample_name}.png', transparent=False, dpi=600)
        logging.info(f"Saved figure to figures/{sample_name}.png.")

    if showMap:
        plt.show()
    else:
        plt.close()
