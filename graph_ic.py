
# first we take a method to read in XYZ files
from modules import arguments
from modules import ic_gen
from modules import out
import networkx as nx
import math



import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



class Atom:
    def __init__(self, symbol, coordinates):
        """
        Constructor for the atom class
        """
        self.symbol = symbol
        self.coordinates = coordinates

    def __repr__(self):
        return f"{self.symbol} at {self.coordinates}"


def readxyz(filename):
    """ 
    Reads in a .xyz file in standard format, returning xyz coordinates and a list of atom_names
    """

    
    xyzf = open(filename, 'r')
    xyzarr = np.zeros([1, 3])
    atomnames = []
    if not xyzf.closed:
        # Read the first line to get the number of particles
        npart = int(xyzf.readline())
        # and next for title card
        title = xyzf.readline()

        # Make an N x 3 matrix of coordinates
        xyzarr = np.zeros([npart, 3])
        i = 0
        for line in xyzf:
            words = line.split()
            if (len(words) > 3):
                atomnames.append(words[0])
                xyzarr[i][0] = float(words[1])
                xyzarr[i][1] = float(words[2])
                xyzarr[i][2] = float(words[3])
                i = i + 1
        molecule = []
        for i in range(npart):
            atom = Atom(atomnames[i],tuple(xyzarr[i]))
            molecule.append(atom)
    return molecule

def list_of_atom_symbols(molecule) -> list:
    """
    Returns a list of just the symbols in the given Molecule

        e.q ['H1','O2','H2']
    """
    symbol_list = []
    for atom in molecule:
        symbol_list.append(atom.symbol)
    return symbol_list 

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

def extract_submolecules(molecule, submolecule_list):
    """ 
    Extracts submolecules from the given molecule

    Attributes:
        molecule: list
            a list of atoms
        submolecule_list: list of list
            a list of list with atom symbols
    """
    submolecules = []

    for submol_symbols in submolecule_list:
        submolecule = [atom for atom in molecule if atom.symbol in submol_symbols]
        submolecules.append(submolecule)
    return submolecules

def draw_graphs(graphs, molecules, figname):
    """
    Uses Networkx to make side-by-side 3D Graphs of multiple molecules.
    `graphs`: list of NetworkX graphs
    `molecules`: list of molecules with corresponding atomic coordinates
    """
    fig = plt.figure(figsize=(len(graphs)*6, 6))  # Make the figure larger to accommodate multiple graphs

    degree_colors = {1: "red", 2: "green", 3: "orange", 4: "lightblue"}

    for idx, (graph, molecule) in enumerate(zip(graphs, molecules), 1):
        ax = fig.add_subplot(1, len(graphs), idx, projection="3d")
        
        # Extract the positions out of the molecule
        pos = {}
        for atom in molecule:
            pos[atom.symbol] = atom.coordinates
        
        # Draw the edges
        for edge in graph.edges():
            x_vals = [pos[edge[0]][0], pos[edge[1]][0]]
            y_vals = [pos[edge[0]][1], pos[edge[1]][1]]
            z_vals = [pos[edge[0]][2], pos[edge[1]][2]]
            ax.plot(x_vals, y_vals, z_vals, color="black")
        
        # Draw the nodes
        degrees = dict(graph.degree())
        x_vals = [pos[node][0] for node in graph.nodes()]
        y_vals = [pos[node][1] for node in graph.nodes()]
        z_vals = [pos[node][2] for node in graph.nodes()]
        node_colors = [degree_colors.get(degrees[node], "green") for node in graph.nodes()]
        ax.scatter(x_vals, y_vals, z_vals, color=node_colors, s=100)
        
        # Label the nodes
        for node, (x, y, z) in pos.items():
            ax.text(x, y, z, f"{node}", size=12, zorder=1, color="k")
        
        # Set labels for axes
        ax.set_xlabel("X axis")
        ax.set_ylabel("Y axis")
        ax.set_zlabel("Z axis")
    
    plt.tight_layout()
    plt.savefig(f"{figname}.png", dpi=600)
    
    #Do show if you want nice animation
    #plt.show()

def asteroidal_triple(graph):
    """ 
    A asteroidal triple is a triple of vertices where two are joined by a path that avoids the neighbors of the third
    """
    print("Hello")

def chain_decomposition(graph):
    """ 
    Peformes the chain decomposition algorithm

    This https://doi.org/10.1016/j.ipl.2013.01.016 tells you all you ever wanted to know about it
    """

    chains = nx.chain_decomposition(graph)
    

def kerighan_lin_algorithm(graph):
    """ 
    Another heuristic algorithm to find partitions

    TODO: This could also be useful but needs work

    This algorithm minimized the number of crossing edges
    """
    print(nx.community.kernighan_lin_bisection(graph))


def connectivity_analysis(graph):
    """  
    This function performes connectivity and cut algorithms on the given molecular graph
    """

    print(nx.junction_tree(graph))

def general_nonlinear(graph,molecule,cov_bonds):
    """
    Set of Functions for Decius Selection
    """
    n_r = len(cov_bonds)
    
    # for n_phi calculate number of nodes with 1 bond
    
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

def general_nonlinear(molecule,graph,cov_bonds):
    """ 
    Calculates the Decius Coordinates for the general nonlinear system
    """

    a = len(molecule)
    b = len(cov_bonds)
    IDOF = len(molecule) * 3 - 6


    n_r = b # num of bonds

    # filter vertices with degree 1
    a_1 = len([node for node,degree in graph.degree() if degree==1])

    n_phi = 4*b - 3*a + a_1
    n_tau = b-a_1 

    return n_r, n_phi, n_tau
    

def ic_generator(graph):
    """ 
    Calculates all the possible IC sets for a given graph
    """

    bonds = el_symm_combinations(find_all_paths_length_n(graph,1))
    angles = el_symm_combinations(find_all_paths_length_n(graph,2))
    dihedrals = el_symm_combinations(find_all_paths_length_n(graph,3))

    return bonds,angles,dihedrals    

def binomial_coefficient(n,k):
    """ 
    Calculates the Binomial Coefficient
    """

    return math.comb(n,k)

    
def main():
    
    args = arguments.get_args()

    if args.xyz:

        # First of read in xyz file

        molecule = readxyz(args.xyz)

    # Initialize out file

    writer = out.Output_Writer("output.txt")
  
    #write header
    writer.write_header(writer.header)

    # First of generate the degree of covalence
    degofc = ic_gen.degree_of_covalance(molecule)


    cov_bonds = ic_gen.covalent_bonds(molecule,degofc)
 
    # With the bonds in place now establish a graph
    G = nx.Graph()
    G.add_nodes_from(list_of_atom_symbols(molecule))
    G.add_edges_from(cov_bonds)
    

    # Next of we make a small graph visualization function
    draw_graphs([G],[molecule], "initial_molecule_graph")

    # Degree Analysis, Calculate Average Parameters
    # Histogramm of degree population
    degree_dict, average_degree,average_degree_neighbour = degree_eval(G)
    
    

    # Matrix Analysis

    _, _, lap_matrix, _ =matrix_analysis(G)


    # Spectral Partitioning --> returns subgraphs

    dis_subgraphs = spectral_bisection(G)

    subgraphs_nodes = []
    for graph in dis_subgraphs:
        subgraphs_nodes.append(list(graph))
    
    # extract submolecules
    submolecules = extract_submolecules(molecule, subgraphs_nodes)

    # Make a picture of the result of spectral partitioning
    draw_graphs(dis_subgraphs,submolecules, "spectral_bisection")


    



    # TODO i need to understand asteroidal structures in the graphs, seems useful
    # Check for asteroidal graph
    if not nx.is_at_free(G):
        print("Hello")

    connectivity_analysis(G)


    


    '''
    Generate possible IC sets
    '''

    bonds,angles,dihedrals = ic_generator(G)

    '''
    Calculate initial combinatorics
    '''
    
    IDOF = 3*len(molecule) -6

    total_combinations = binomial_coefficient(len(bonds)+len(angles)+len(dihedrals),IDOF)


    print(len(angles))
    writer.write_subheader("== Total Combinations ICs ===")
    writer.write_line(f"(IDOF, TOT_ICs) = {total_combinations}")

    

     

    # Decius Calculations
    
 


    




if __name__ == "__main__":
    main()