# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import scipy.stats
import pandas as pd
from jmct import time_elapsed
from functools import cache

@cache
def findClusterProbability(
        cluster_length: int, 
        number_of_mutations: int, 
        per_bp_mut_probability: float
        ) -> float:
    """
    Calculates the probability of observing a cluster of mutations of a given length 
    and with a given number of mutations, given the probability of a mutation occurring at any given base pair.

    Args:
        cluster_length (int): The length of the cluster to calculate the probability for.
        number_of_mutations (int): The total number of mutations observed in the cluster.
        per_bp_mut_probability (float): The probability of a mutation occurring at any given base pair.

    Returns:
        float: The probability of observing a cluster of mutations of a given length
        and with a given number of mutations, given the probability of a mutation occurring at any given base pair.
    """
    cdf = scipy.stats.nbinom.cdf(cluster_length, number_of_mutations, per_bp_mut_probability)
    return cdf

def computeClusterSignificanceOfAllMutationPairs(
        possible_cluster_mutation_list: list, 
        per_site_probability: float
        ) -> list:
    """
    Computes the significance of all possible cluster lengths for a given list of mutations.

    Args:
        possible_cluster_mutation_list (list): A list of possible mutations to consider.
        per_site_probability (float): The probability of a mutation occurring at any given base pair.

    Returns:
        list: A list of tuples containing the start and end positions of each cluster, the number of mutations in the cluster, and the probability of observing a cluster of this length with this many mutations.
    """

    #compute all pairs of values in possible_cluster_mutation_list
    possible_pairs = [ (possible_cluster_mutation_list[i], possible_cluster_mutation_list[j]) for i in range(len(possible_cluster_mutation_list)) for j in range(i+1, len(possible_cluster_mutation_list)) ]

    cluster_probabilities = []
    for pair in possible_pairs:
        cluster_length = pair[1] - pair[0]
        
        #compute the number of elements in possible_cluster_mutation_list that are in the range [pair[0], pair[1]]
        number_of_mutations = len([x for x in possible_cluster_mutation_list if x >= pair[0] and x <= pair[1]])

        p_cdf = findClusterProbability(
            cluster_length = cluster_length, 
            number_of_mutations = number_of_mutations, 
            per_bp_mut_probability = per_site_probability)
        
        #cluster_probabilities contains a tuple with:
            # the start of the cluster
            # the end of the cluster
            # the number of mutations in the cluster
            # the probability of observing a cluster of this length with this many mutations
        cluster_probabilities.append((pair[0], pair[1], number_of_mutations, p_cdf))

    return cluster_probabilities

def findPotentialClusteredMutations(
        mutation_list: list, 
        max_distance_between_mutations: int = 10000
        ) -> list[list]:
    """
    Splits a list of mutations into sublists where the distance between two mutations is less than or equal to a given maximum distance.
    Sublists with less than 2 elements are discarded since 1 mutation cannot form a cluster.

    Args:
        mutation_list (list): A list of mutations to split into sublists.
        max_distance_between_mutations (int): The maximum distance between two mutations to include them in the same sublist. Defaults to 10000.

    Returns:
        list[list]: A list of sublists, where each sublist contains mutations that are within the maximum distance of each other.
        Sublists with less than 2 elements are discarded. 
    """
    potential_clusters = []
    #go through the mutation list and split it into sublists at the points where the distance between two mutations is greater than max_distance_between_mutations
    #start with the second mutation in the list and calculate the distance between it and the previous mutation
    #if the distance is smaller than max_distance_between_mutations, add both mutations to the current sublist and continue to the next mutation
    #if the distance is greater than max_distance_between_mutations, close the current sublist and start a new one, then continue to the next mutation

    current_sublist = [mutation_list[0]]

    # Below replaced with further for efficiency
    # for i in range(1, len(mutation_list)):

    #     if mutation_list[i] - mutation_list[i-1] <= max_distance_between_mutations:
    #         current_sublist.append(mutation_list[i])
    #     else:
    #         potential_clusters.append(current_sublist)
    #         current_sublist = []
    #         current_sublist.append(mutation_list[i])
    
    # potential_clusters.append(current_sublist)

    # #go through the list of sublists and remove all sublists that have less than 2 elements
    # potential_clusters = [x for x in potential_clusters if len(x) >= 2]
    
    for i in range(1, len(mutation_list)):
        if mutation_list[i] - mutation_list[i-1] <= max_distance_between_mutations:
            current_sublist.append(mutation_list[i])
        else:
            if len(current_sublist) >= 2:
                potential_clusters.append(current_sublist)
            current_sublist = [mutation_list[i]]
    
    if len(current_sublist) >= 2:
        potential_clusters.append(current_sublist)


    return potential_clusters

def find_largest_significant_clusters(
        potential_clusters: list[list], 
        per_site_probability: float, 
        min_p_value: float
        ) -> list:
    """
    Finds the largest clusters of mutations in a list of potential clusters that have a p-value less than a given threshold.

    Args:
        potential_clusters (list[list]): A list of sublists, where each sublist contains mutations that are within a certain minimum distance of each other.
        per_site_probability (float): The probability of a mutation occurring at any given base pair.
        min_p_value (float): The minimum p-value required for a cluster to be considered significant.

    Returns:
        list[tuple]: A list of tuples containing the start and end positions of each significant cluster, the number of mutations in the cluster, and the probability of observing a cluster of this length with this many mutations.
    """
    significant_clusters = []
    for possible_cluster_mutation_list in potential_clusters:
        cluster_probabilities = computeClusterSignificanceOfAllMutationPairs(
            possible_cluster_mutation_list, 
            per_site_probability)
        # find the longest cluster with a p-value less than min_p_value
        longest_cluster = None
        for cluster in cluster_probabilities:
            if cluster[3] < min_p_value:
                if longest_cluster is None:
                    longest_cluster = cluster
                elif cluster[1] - cluster[0] > longest_cluster[1] - longest_cluster[0]:
                    longest_cluster = cluster

        if longest_cluster is not None:
            significant_clusters.append(longest_cluster)

    return significant_clusters

def reportSignificantClusters(
        significant_clusters: list[tuple], 
        min_p_value: float, 
        per_site_probability: float
        ) -> None:
    """
    Prints out a report of the significant clusters of mutations found. Debugging function.

    Args:
        significant_clusters (list[tuple]): A list of tuples containing the start and end positions of each significant cluster, the number of mutations in the cluster, and the probability of observing a cluster of this length with this many mutations.
        min_p_value (float): The minimum p-value required for a cluster to be considered significant.
        per_site_probability (float): The probability of a mutation occurring at any given base pair.

    Returns:
        None: This function does not return anything, it simply prints the report to the console.
    """
    print(f"Significant clusters (p<{min_p_value}) with per site mutation probability {per_site_probability:.6e}:")
    for significant_cluster in significant_clusters:
        print(f"    {significant_cluster[0]:>8} - {significant_cluster[1]:<8} ({(significant_cluster[1]-significant_cluster[0]):>5} bp) with {significant_cluster[2]:3} mutations (p={significant_cluster[3]:.6e})")

def findSampleClusters(
        df_sample: pd.DataFrame, 
        per_site_probability: float, 
        min_p_value: float, 
        max_distance_between_mutations: int = 10000,
        report=False) -> pd.DataFrame:
    """
    Finds clusters of mutations in a given DataFrame of mutation data.

    Args:
        df_sample (pd.DataFrame): A DataFrame containing mutation data, with columns 'Chromosome' and 'Region'.
        per_site_probability (float): The probability of a mutation occurring at any given base pair.
        min_p_value (float): The minimum p-value required for a cluster to be considered significant.
        max_distance_between_mutations (int): The maximum distance between two mutations to include them in the same cluster. Defaults to 10000.
        report (bool): Whether or not to print a report of the significant clusters found. Defaults to False.

    Returns:
        pd.DataFrame: A DataFrame containing information about the significant clusters found, with columns 'Chromosome', 'Start', 'End', 'Length', 'Mutations', 'P_value', 'Per_site_mutation_probability', and 'Cluster_P_value_Threshold'.
    """
    #get the list of chromosomes
    sample_chromosomes = df_sample['Chromosome'].unique()

    df_clusters = pd.DataFrame(
        columns=[
            'Chromosome', 
            'Start', 
            'End', 
            'Length', 
            'Mutations', 
            'P_value', 
            'Per_site_mutation_probability', 
            'Cluster_P_value_Threshold'])
    
    #go through each chromosome and find potential clusters of mutations
    for sample_chromosome in sample_chromosomes:
        if sample_chromosome == 'ref|NC_001224|':
            continue
        #get the list of mutations positions for the current chromosome, as integers
        df_sample['Region'] = df_sample['Region'].astype(int)
        mutation_list = df_sample[df_sample['Chromosome'] == sample_chromosome]['Region'].tolist()

        #sort the list of mutations
        mutation_list.sort()

        #find potential clusters of mutations
        potential_clusters = findPotentialClusteredMutations(
            mutation_list, 
            max_distance_between_mutations = max_distance_between_mutations)
        
        significant_clusters = find_largest_significant_clusters(
            potential_clusters, 
            per_site_probability, 
            min_p_value)

        if report and significant_clusters != []:
            print(f"Chromosome {sample_chromosome}")
            reportSignificantClusters(
                significant_clusters, 
                min_p_value, 
                per_site_probability)

        #add the significant clusters to the dataframe
        for significant_cluster in significant_clusters:
            df_clusters = pd.concat([
                df_clusters,
                pd.DataFrame(
                    [[
                        sample_chromosome, 
                        significant_cluster[0], 
                        significant_cluster[1], 
                        significant_cluster[1]-significant_cluster[0], 
                        significant_cluster[2], 
                        significant_cluster[3], 
                        per_site_probability, 
                        min_p_value]], 
                    columns=[
                        'Chromosome',
                        'Start', 
                        'End', 
                        'Length', 
                        'Mutations', 
                        'P_value', 
                        'Per_site_mutation_probability', 
                        'Cluster_P_value_Threshold']
                    )],
                
                ignore_index=True)

    return(df_clusters)

@time_elapsed
def analyzeSampleForMutationClusters(
        df_mutations: pd.DataFrame, 
        genome_size: int,
        min_p_value: float = 0.001,
        max_distance_between_mutations: int = 10000,
        min_cluster_mutations: int = 2
        ) -> pd.DataFrame:
    """
    Analyzes a given DataFrame of mutation data for significant clusters of mutations.

    Args:
        df_mutations (pd.DataFrame): A DataFrame containing mutation data, with columns 'Chromosome' and 'Region'.
        genome_size (int): The size of the genome being analyzed (in base pairs).
        min_p_value (float): The minimum p-value required for a cluster to be considered significant. Defaults to 0.001.
        max_distance_between_mutations (int): The maximum distance between two mutations to include them in the same cluster. Defaults to 10000.
        min_cluster_mutations (int): The minimum number of mutations required for a cluster to be considered significant. Defaults to 2.

    Returns:
        pd.DataFrame: A DataFrame containing information about the significant clusters found, with columns 'Chromosome', 'Start', 'End', 'Length', 'Mutations', 'P_value', 'Per_site_mutation_probability', and 'Cluster_P_value_Threshold'.
    """
    total_outcome_mutations = len(df_mutations)
    per_site_probability = total_outcome_mutations/genome_size

    df_clusters = findSampleClusters(
        df_mutations, 
        per_site_probability, 
        min_p_value, 
        max_distance_between_mutations)
    
    df_clusters = df_clusters[df_clusters['Mutations'] >= min_cluster_mutations]
    df_clusters.reset_index(drop=True, inplace=True)

    return df_clusters
