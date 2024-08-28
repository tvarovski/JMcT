# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

# Go through the entire file and make sure that all settings are correct for your data

settings = {
    'data_path': "data",
    'pval_list': [0.01, 0.001, 0.0001, 0.00001],
    'max_clust_distance': 10000,
    'min_clust_mutations': 3,
    'sample_name_len': 8 # number sample file leading charachters used for output names. Trailing characters will be omitted. Set to 0 to use entire name
    }

# The following list of chromosomes is based on the S288C yeast reference genome
# The list is ordered by chromosome name and the length of the chromosome in base pairs
# make sure to update the list if you are using a different reference genome
chromosomes: list[list] = [
    ['ref|NC_001133|',230218], #1.91% -> % of the entire genome, not essential (hence in comments)
    ['ref|NC_001134|',813184], #6.74%
    ['ref|NC_001135|',319343], #2.64%
    ['ref|NC_001136|',1531933],#12.69%
    ['ref|NC_001137|',576070], #4.77%
    ['ref|NC_001138|',270161], #2.24%
    ['ref|NC_001139|',1090940],#9.03%
    ['ref|NC_001140|',562643], #4.66%
    ['ref|NC_001141|',439888], #3.64%
    ['ref|NC_001142|',745751], #6.17%
    ['ref|NC_001143|',666816], #5.52%
    ['ref|NC_001144|',1078177],#8.93%
    ['ref|NC_001145|',924431], #7.65%
    ['ref|NC_001146|',784333], #6.49%
    ['ref|NC_001147|',1091291],#9.03%
    ['ref|NC_001148|',948066]] #7.85%


# This is a list of centromere locations for the S288C yeast reference genome
# chromosome names need to be numeric, not roman numerals or any other format
# it is needed for figure generation
import pandas as pd
from numpy import NaN
chr_centromeres = pd.DataFrame([
    ["chr1" , 151582, "cen", NaN, NaN, NaN, "cen"],
    ["chr2" , 238323, "cen", NaN, NaN, NaN, "cen"],
    ["chr3" , 114501, "cen", NaN, NaN, NaN, "cen"],
    ["chr4" , 449821, "cen", NaN, NaN, NaN, "cen"],
    ["chr5" , 152104, "cen", NaN, NaN, NaN, "cen"],
    ["chr6" , 148627, "cen", NaN, NaN, NaN, "cen"],
    ["chr7" , 497038, "cen", NaN, NaN, NaN, "cen"],
    ["chr8" , 105703, "cen", NaN, NaN, NaN, "cen"],
    ["chr9" , 355745, "cen", NaN, NaN, NaN, "cen"],
    ["chr10", 436425, "cen", NaN, NaN, NaN, "cen"],
    ["chr11", 440246, "cen", NaN, NaN, NaN, "cen"],
    ["chr12", 150947, "cen", NaN, NaN, NaN, "cen"],
    ["chr13", 268149, "cen", NaN, NaN, NaN, "cen"],
    ["chr14", 628875, "cen", NaN, NaN, NaN, "cen"],
    ["chr15", 326702, "cen", NaN, NaN, NaN, "cen"],
    ["chr16", 556073, "cen", NaN, NaN, NaN, "cen"]], 

    columns=["Chromosome", "Region", "Reference", "Allele", "Frequency", "Zygosity", "SPECTRA_STRANDWISE"])

#### DO NOT MODIFY BELOW THIS LINE unless you know what you are doing ####

am1003_to_roman : dict[str,str] = {
    'chr_1_AM_1003_v1':"I",
    'chr_2_AM_1003_v1':"II",
    'chr_3_AM_1003_v1':"III",
    'chr_4_AM_1003_v1':"IV",
    'chr_5_AM_1003_v1':"V",
    'chr_6_AM_1003_v1':"VI",
    'chr_7_AM_1003_v1':"VII",
    'chr_8_AM_1003_v1':"VIII",
    'chr_9_AM_1003_v1':"IX",
    'chr_10_AM_1003_v1':"X",
    'chr_11_AM_1003_v1':"XI",
    'chr_12_AM_1003_v1':"XII",
    'chr_13_AM_1003_v1':"XIII",
    'chr_14_AM_1003_v1':"XIV",
    'chr_15_AM_1003_v1':"XV",
    'chr_16_AM_1003_v1':"XVI"
    }

am1003_to_numeric : dict[str,str] = {
    'chr_1_AM_1003_v1':"chr1",
    'chr_2_AM_1003_v1':"chr2",
    'chr_3_AM_1003_v1':"chr3",
    'chr_4_AM_1003_v1':"chr4",
    'chr_5_AM_1003_v1':"chr5",
    'chr_6_AM_1003_v1':"chr6",
    'chr_7_AM_1003_v1':"chr7",
    'chr_8_AM_1003_v1':"chr8",
    'chr_9_AM_1003_v1':"chr9",
    'chr_10_AM_1003_v1':"chr10",
    'chr_11_AM_1003_v1':"chr11",
    'chr_12_AM_1003_v1':"chr12",
    'chr_13_AM_1003_v1':"chr13",
    'chr_14_AM_1003_v1':"chr14",
    'chr_15_AM_1003_v1':"chr15",
    'chr_16_AM_1003_v1':"chr16"
    }

roman_to_numeric : dict[str,str] = {
    "I": "chr1",
    "II": "chr2",
    "III": "chr3",
    "IV": "chr4",
    "V": "chr5",
    "VI": "chr6",
    "VII": "chr7",
    "VIII": "chr8",
    "IX": "chr9",
    "X": "chr10",
    "XI": "chr11",
    "XII": "chr12",
    "XIII": "chr13",
    "XIV": "chr14",
    "XV": "chr15",
    "XVI": "chr16"
    }

numeric_to_roman : dict[str,str] = {
    "chr1" : "I",
    "chr2" : "II",
    "chr3" : "III",
    "chr4" : "IV",
    "chr5" : "V",
    "chr6" : "VI",
    "chr7" : "VII",
    "chr8" : "VIII",
    "chr9" : "IX",
    "chr10": "X",
    "chr11": "XI",
    "chr12": "XII",
    "chr13": "XIII",
    "chr14": "XIV",
    "chr15": "XV",
    "chr16": "XVI"
    }

S288C_to_numeric : dict[str,str] = {
    'ref|NC_001133|': "chr1",
    'ref|NC_001134|': "chr2",
    'ref|NC_001135|': "chr3",
    'ref|NC_001136|': "chr4",
    'ref|NC_001137|': "chr5",
    'ref|NC_001138|': "chr6",
    'ref|NC_001139|': "chr7",
    'ref|NC_001140|': "chr8",
    'ref|NC_001141|': "chr9",
    'ref|NC_001142|': "chr10",
    'ref|NC_001143|': "chr11",
    'ref|NC_001144|': "chr12",
    'ref|NC_001145|': "chr12",
    'ref|NC_001146|': "chr14",
    'ref|NC_001147|': "chr15",
    'ref|NC_001148|': "chr16"
    }

S288C_to_roman : dict[str,str] = {
    'ref|NC_001133|': "I",
    'ref|NC_001134|': "II",
    'ref|NC_001135|': "III",
    'ref|NC_001136|': "IV",
    'ref|NC_001137|': "V",
    'ref|NC_001138|': "VI",
    'ref|NC_001139|': "VII",
    'ref|NC_001140|': "VIII",
    'ref|NC_001141|': "IX",
    'ref|NC_001142|': "X",
    'ref|NC_001143|': "XI",
    'ref|NC_001144|': "XII",
    'ref|NC_001145|': "XIII",
    'ref|NC_001146|': "XIV",
    'ref|NC_001147|': "XV",
    'ref|NC_001148|': "XVI"
    }
