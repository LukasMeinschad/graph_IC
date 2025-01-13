####################################################################################################
# Description: This module contains functions for graph analysis.
####################################################################################################

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def matrix_analysis(graph):
    """ 
    Makes adjecency, incidence and laplacian matrix 


    Proposititon (because i will sure forget it):

    Let G be an undirected graph, the multiplicity of k of the eigenvalue 0 of 
    the laplacian matrix equals the number of connected components

    
    Returns:
        adj_matrix, lap_matrix, spectrum_lap, connected_components
    """
    adj_matrix = nx.adjacency_matrix(graph)
    lap_matrix = nx.laplacian_matrix(graph)
    spectrum_lap = nx.laplacian_spectrum(graph)
    
    # Set values < 10⁻16 to zero

    spectrum_lap[spectrum_lap < 10E-16] = 0
    
    # Determine the connected components with the zero eigenvalues
    connected_components = np.count_nonzero(spectrum_lap==0)

    return adj_matrix, lap_matrix, spectrum_lap, connected_components

def degree_eval(graph):
    """ 
    Evaluates the degree and stores them to a dictionary and calculates the
    average degree of a graph

    Attributes:
        graph: A networkx graph
    
    Returns:
        degree_dict: A dictionary with all degrees of the nodes
        average_degree: The average degree of the graph
        average_degree_neighbour:
    """

    degree_dict = {}

    for node in graph.nodes():
        degree_dict[f"{node}"] = graph.degree[node]
    
    average_degree_neighbour = nx.average_neighbor_degree(graph)

    # Calculate the average degree

    # Number of nodes
    abs_V = graph.number_of_nodes()

    sum_degree = 0
    for node in graph.nodes():
        sum_degree += graph.degree[node]
    average_degree = sum_degree/abs_V

    # Further we generate a histogramm of the degree population

    degree_sequence = sorted((d for n, d in graph.degree()),reverse=True)
    
    plt.figure(figsize=(6,6))
    plt.bar(*np.unique(degree_sequence, return_counts = True))
    plt.title("Degree Histogramm")
    plt.xlabel("Degree")
    plt.ylabel("# of Nodes")
    plt.savefig("degree_histogramm.png",dpi=600)

    return degree_dict, average_degree, average_degree_neighbour

def dominating_set(G):
    """
    This function calculates the dominating set of a graph G

    A dominating set of a graph G is a subset D of V such that every vertex not in D is adjacent to at least one member of D

    Attributes:
        G: A networkx graph
    
    Returns:
        dominating_set: A list with the dominating set
    """

    dominating_set = nx.dominating_set(G)

    subgraph = G.subgraph(dominating_set)

    return list(dominating_set), subgraph