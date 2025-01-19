############################################
# This module containts the distrubution of ICs to subgraphs, fixation algorithms and subset combinations

import numpy as np
import networkx as nx

from . import decius_rules 
from . import general_functions

def spectral_bisection(graph):
    """ 
    This method uses the Fielder Vector to bisect the graph
    
    So the Fielder uses the second eigenvector of the Laplacian matrix and can be used to bisect a graph

    The algorithm is as follows:

        1. Compute Fielder Eigenvector
        2. Search for median of v
        for each node check if v[i] ≤ median
            put node in partition V1
        else
            put node in partition V2

    """
    dissected_graph = nx.spectral_bisection(graph)
    
    list_of_subgraphs = []
    for tup in dissected_graph:
        list_of_subgraphs.append(graph.subgraph(tup))

    return list_of_subgraphs


def dfs_path(G,start_node, path_length):
    """ 
    Performes a DFS algorithm searching for paths with a given length
    """
    def dfs(current_node, current_path):
        # If the current path has the required length (path_length + 1)
        # Append it to the possible paths
        if len(current_path) == path_length + 1:
            paths.append(list(current_path))
        

        # next recursively visit all the neighbours
        for neighbor in G.neighbors(current_node):
            if neighbor not in current_path: # avoid the revisiting
                current_path.append(neighbor)
                dfs(neighbor,current_path)
                current_path.pop() # backtrack
    
    paths = []
    dfs(start_node,[start_node])
    return paths

def find_all_paths_length_n(G,n):
    """ 
    This function then find all paths
    """
    all_paths = []
    for node in G.nodes():
        # we iterate over all the nodes
        node_paths = dfs_path(G,node,n)
        all_paths.extend(node_paths)
    return all_paths


def el_symm_combinations(input_list):
    # Sort
    seen = set()
    unique_combinations = []

    for sublist in input_list:
        sorted_tuple = (tuple(sorted(sublist)))
        if sorted_tuple not in seen:
            seen.add(sorted_tuple)
            unique_combinations.append(sublist)
    return unique_combinations



def ic_generator(graph):
    """ 
    Calculates all the possible IC sets for a given graph
    """

    bonds = el_symm_combinations(find_all_paths_length_n(graph,1))
    angles = el_symm_combinations(find_all_paths_length_n(graph,2))
    dihedrals = el_symm_combinations(find_all_paths_length_n(graph,3))

    return bonds,angles,dihedrals    



def subgraph_ic_fixation(subgraph,submolecule):

    # Generate all ICs of the Subgraph
    bonds_sub, angles_sub, dihedrals_sub = ic_generator(subgraph)

    bonds_sub = [tuple(elem) for elem in bonds_sub]
    angles_sub = [tuple(elem) for elem in angles_sub]
    dihedral_sub = [tuple(elem) for elem in dihedrals_sub]
 
    
    n_r_sub, n_phi_sub, n_tau_sub = decius_rules.general_nonlinear(submolecule,subgraph,bonds_sub)

    # TODO WRITER for out


    # Calculate the number of combinations 
    bond_combs_sub = general_functions.binomial_coefficient(len(bonds_sub),n_r_sub)
    if n_phi_sub > 0:
        angle_combs_sub = general_functions.binomial_coefficient(len(angles_sub),n_phi_sub)
    else:
        angle_combs_sub = 1
    if n_tau_sub > 0: 
        dihedral_combs_sub = general_functions.binomial_coefficient(len(dihedrals_sub),n_tau_sub)
    else:
        dihedral_combs_sub = 1



    # Initialize the Nomodeco IC dict

    ic_dict_temp = {
        "bonds": [],
        "angles": [],  
        "linear valence angles": [],
        "out of plane angles": [],
        "dihedrals": [], 
    } 

    # initialize a empty dictionary with the total combinations

    


