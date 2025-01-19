####################################################################################################
# Description: This module contains functions for graph analysis.
####################################################################################################

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from scipy.cluster.hierarchy import dendrogram, linkage


# This is used for Kernighan_lin_Bisection
from networkx.algorithms import community


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

def wiener_index(G):
    """
    Calculates the wiener_index for a given graph
    """
    return nx.wiener_index(G)

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
    dominating_set = subgraph.nodes()

    return dominating_set, subgraph

def graph_reduction(G):
    """ 
    Function that removes all vertices with degree 0 or 1 from the graph.
    This is a pruning step for fixing internal coordinates 
    """

    # Copy the graph to not alter original memory allocation
    copy_graph = G.copy()
    while True:
        low_degree_nodes = [node for node in copy_graph.nodes() if G.degree(node) < 2]
        if not low_degree_nodes:
            break
        copy_graph.remove_nodes_from(low_degree_nodes)
    return copy_graph

# Agglomerative Clustering

def agglomerative_clustering(G,n_clusters=3):
    """
    Performs agglomerative clustering on the graph
    Default Clusters are 2
    """
    # a adjacency matrix is needed for the clustering

    adj_matrix = nx.to_numpy_array(G)

    clustering = AgglomerativeClustering(n_clusters=n_clusters).fit(adj_matrix)
    labels = clustering.fit_predict(1-adj_matrix)

    # Compute the linkage matrix
    Z = linkage(adj_matrix, 'ward')

    # Plot the dendrogram
    plt.figure(figsize=(10, 7))
    dendrogram(Z,labels=list(G.nodes()))
    plt.title("Dendrogram for Agglomerative Clustering")
    plt.xlabel("Node Index")
    plt.ylabel("Distance")
    plt.savefig("dendrogram.png", dpi=600)
    #plt.show()

    # Print the Clustering
    pos = nx.spring_layout(G)
    plt.figure(figsize=(6,6))
    nx.draw(G, pos, node_color=labels, with_labels=True)
    plt.title("Agglomerative Clustering")
    plt.savefig("agglomerative_clustering.png",dpi=600)

    # maybe not the best approach we need connected graphs



def k_means_clustering(G,n_clusters=2):
    """
    Performs k-means clustering on the graph
    Default Clusters are 2
    """

    # a adjacency matrix is needed for the clustering

    adj_matrix = nx.to_numpy_array(G)

    clustering = KMeans(n_clusters=n_clusters).fit(adj_matrix)
    labels = clustering.fit_predict(1-adj_matrix)

    # Print the Clustering
    pos = nx.spring_layout(G)
    plt.figure(figsize=(6,6))
    nx.draw(G, pos, node_color=labels, with_labels=True)
    plt.title("K-Means Clustering")
    plt.savefig("k_means_clustering.png",dpi=600)

    # maybe not the best approach we need connected graphs


def dbscan_clustering(G, eps=0.1, min_samples=2):
    """
    Performs DBSCAN clustering on a graph.

    Args:
        G (networkx.Graph): The input graph.
        eps (float): The maximum distance between two samples for one to be considered as in the neighborhood of the other.
        min_samples (int): The number of samples in a neighborhood for a point to be considered as a core point.

    Returns:
        labels (list): Cluster labels for each node.
    """
    # Convert the graph to an adjacency matrix
    adjacency_matrix = nx.to_numpy_array(G)

    # Perform DBSCAN clustering
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = dbscan.fit_predict(1 - adjacency_matrix)
    pos = nx.spring_layout(G)
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, node_color=labels, with_labels=True, cmap=plt.cm.rainbow, node_size=500, font_size=10)
    plt.title("DBSCAN Clustering of the Graph")
    plt.savefig("dbscan_clustering.png", dpi=600)

    
    return labels


def kerighan_lin_algorithm(G):
    """ 
    Input is undirected Graph G =(V,E) with vertex set V, edge set E

    Graph is partitioned into two disjoint sets A, B
    """
    return community.kernighan_lin_bisection(G) 
