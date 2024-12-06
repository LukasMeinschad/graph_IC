
# first we take a method to read in XYZ files
from modules import arguments
from modules import ic_gen
import networkx as nx


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


def draw_graph(graph,molecule):
    """ 
    Uses Networkx to make a cool 3D Graph of your grah
    """
    # make a dictionary to store the node colors
    node_colors = []
    for node in graph.nodes():
        if graph.degree[node] == 2:
            node_colors.append("red")
        elif graph.degree[node] == 3:
            node_colors.append("lightblue")
        elif graph.degree[node] == 4:
            node_colors.append("orange")
        else:
            node_colors.append("green")
    
    node_labels = {node: f"{node}" for node in graph.nodes()}

    print(graph.nodes())
    print(molecule)
    # Extract the positions out of the molecule
    pos = {}
    i= 0
    for atom in molecule:
        pos[atom.symbol] = atom.coordinates
        i +=1
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # This draws the edges
    for edge in graph.edges():
        x_vals = [pos[edge[0]][0],pos[edge[1]][0]]
        y_vals = [pos[edge[0]][1],pos[edge[1]][1]]
        z_vals = [pos[edge[0]][2],pos[edge[1]][2]]
        ax.plot(x_vals,y_vals,z_vals,color="black")

   # And we hilight the degrees of the individual nodes

    degree_colors = {1: "red", 2: "green", 3: "orange", 4: "lightblue"}

    degrees = dict(graph.degree())

     
    # Now draw the nodes
    x_vals = [pos[node][0] for node in graph.nodes()]
    y_vals = [pos[node][1] for node in graph.nodes()]
    z_vals = [pos[node][2] for node in graph.nodes()]
    
    node_colors = [degree_colors[degrees[node]] for node in graph.nodes()]
    ax.scatter(x_vals,y_vals,z_vals, color =node_colors,s=100)

    # Label the nodes 
    for node, (x,y,z) in pos.items():
        ax.text(x,y,z, f"{node}",size=12, zorder=1, color="k")

    # Label the Degrees
    legend_elements = [
        plt.Line2D([0],[0], marker="o", color = "w", label="Degree 1",markeredgecolor="red"),
        plt.Line2D([0],[0], marker="o", color = "w", label="Degree 2",markeredgecolor="green"),
        plt.Line2D([0],[0], marker="o", color = "w", label="Degree 2",markeredgecolor="orange"),
        plt.Line2D([0],[0], marker="o", color = "w", label="Degree 1",markeredgecolor="lightblue")
    ]

    ax.legend(handles=legend_elements)
        

    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")
    
    
    plt.savefig("initial_molecular_graph.png",dpi=600)

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

    """
    adj_matrix = nx.adjacency_matrix(graph)
    lap_matrix = nx.laplacian_matrix(graph)
    spectrum_lap = nx.laplacian_spectrum(graph)
    
    # Set values < 10⁻16 to zero

    spectrum_lap[spectrum_lap < 10E-16] = 0
    
    # Determine the connected components with the zero eigenvalues
    connected_components = np.count_nonzero(spectrum_lap==0)

    
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

def main():
    
    args = arguments.get_args()

    if args.xyz:

        # First of read in xyz file

        molecule = readxyz(args.xyz)

    # First of generate the degree of covalence
    degofc = ic_gen.degree_of_covalance(molecule)


    cov_bonds = ic_gen.covalent_bonds(molecule,degofc)
 
    # With the bonds in place now establish a graph
    G = nx.Graph()
    G.add_nodes_from(list_of_atom_symbols(molecule))
    G.add_edges_from(cov_bonds)
    

    # Next of we make a small graph visualization function
    draw_graph(G,molecule)

    # Degree Analysis, Calculate Average Parameters
    # Histogramm of degree population
    degree_dict, average_degree,average_degree_neighbour = degree_eval(G)
    
    # Adjacency Matrix

    matrix_analysis(G)


    # First simple algorithm is to find all possible paths

    # TODO eliminate the symmetric tuples here
    
    paths_lenght_2 = find_all_paths_length_n(G,2)
    paths_length_3 = find_all_paths_length_n(G,3)

    


    




if __name__ == "__main__":
    main()