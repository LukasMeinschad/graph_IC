
# first we take a method to read in XYZ files
from modules import arguments
from modules import ic_gen
from modules import out
from modules import graph_visualization
from modules import graph_analysis
from modules import general_functions
from modules import ic_distribution 


import networkx as nx
import math
import itertools


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
    



def connectivity_analysis(graph):
    """  
    This function performes connectivity and cut algorithms on the given molecular graph
    """

    print(nx.junction_tree(graph))


 
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
    
    Attributes:
        molecule:
            a list of atoms
        graph:
            a networkx graph
        cov_bonds:
            a list of covalent bonds
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
    graph_visualization.draw_graphs([G],[molecule], "initial_molecule_graph")

    '''
    General Graph Analysis
    ''' 

    chemical_formula = general_functions.chemical_formula(molecule)
     

    degree_dict, average_degree,average_degree_neighbour = graph_analysis.degree_eval(G)
    adj_matrix, lap_matrix, spectrum_lap, connected_components = graph_analysis.matrix_analysis(G)
    wiener_index = graph_analysis.wiener_index(G)



    dominating_set, dominating_set_graph= graph_analysis.dominating_set(G)

    dominating_set_atoms = (extract_submolecules(molecule,dominating_set))
    
    # Weird Fix this some time
    flattened_dominating_set_atoms = [atom for sublist in dominating_set_atoms for atom in sublist]

    graph_visualization.draw_graphs([dominating_set_graph],[flattened_dominating_set_atoms], "dominating_set")

    writer.write_subheader("=== General Graph Analysis ===")
    writer.write_line("Chemical Formula: " + chemical_formula)
    writer.write_line("Average Degree: " + f"{average_degree:.2f}")
    writer.write_line("Wiener Index: " + f"{wiener_index}")


    

    
    '''
    Clustering Approaches
    '''

    # Agglomerative Clustering

    
    #graph_analysis.agglomerative_clustering(G)
    #graph_analysis.k_means_clustering(G)
    #graph_analysis.dbscan_clustering(G)

    # Graph Pruning:
    # Ref : Graphical Approach for Defining Natural Internal Coordinates
    # Basically remove all nodes with degree 1 or zero
    pruned_subgraph = graph_analysis.graph_reduction(G)
    # now we need to extract the submolecules
    pruned_submolecules = extract_submolecules(molecule,[list(pruned_subgraph)])
    # flatten the submolecules
    pruned_submolecules = [atom for sublist in pruned_submolecules for atom in sublist]
    graph_visualization.draw_graphs([pruned_subgraph],[pruned_submolecules], "pruned_graph")



    dis_subgraphs = spectral_bisection(G)

   
    subgraphs_nodes = []
    for graph in dis_subgraphs:
        subgraphs_nodes.append(list(graph))

    
    # extract submolecules
    submolecules = extract_submolecules(molecule, subgraphs_nodes)


    # Make a picture of the result of spectral partitioning
    graph_visualization.draw_graphs(dis_subgraphs,submolecules, "spectral_bisection")
    
    # TODO i need to understand asteroidal structures in the graphs, seems useful
    # Check for asteroidal graph
    if not nx.is_at_free(G):
        print("Hello")

    connectivity_analysis(G)

    '''
    Generate possible IC sets
    '''

    bonds,angles,dihedrals = ic_generator(G)

    # TODO fix this

    bonds = [tuple(elem) for elem in bonds]
    angles= [tuple(elem) for elem in angles]
    dihedrals = [tuple(elem) for elem in dihedrals]
    writer.write_subheader("=== ICs generated for Molecule ===")
    writer.write_line("Bonds")
    writer.write_data(bonds)
    writer.write_line("Angles")
    writer.write_data(angles)
    writer.write_line("Dihedrals")
    writer.write_data(dihedrals)
   

    '''
    Calculate initial combinatorics
    '''
    
    IDOF = 3*len(molecule) -6

    total_combinations = binomial_coefficient(len(bonds)+len(angles)+len(dihedrals),IDOF)

      
    writer.write_subheader("=== Total Combinations ICs ===")
    writer.write_line(f"(IDOF, TOT_ICs) = {total_combinations}")

    '''
    Apply Decius Rules to total molecule
    '''

    n_r, n_phi, n_tau = general_nonlinear(molecule,G,cov_bonds)

    
    writer.write_subheader("=== IC Types needed (DECIUS) ===")
    writer.write_line(f"Bonds needed: {n_r}")
    writer.write_line(f"Angles needed: {n_phi}")
    writer.write_line(f"Dihedrals needed: {n_tau}")

    # Calculate the possible combinations 

    bond_combs = binomial_coefficient(len(bonds), n_r)
    angle_combs = binomial_coefficient(len(angles),n_phi)
    dihedral_combs = binomial_coefficient(len(dihedrals),n_tau)

    writer.write_line(f"Bond subsets: {bond_combs}")
    writer.write_line(f"Angle subsets: {angle_combs}")
    writer.write_line(f"Dihderal subsets: {dihedral_combs}")
    writer.write_line("\n")

    total_ic_sets = bond_combs * angle_combs * dihedral_combs

    writer.write_line(f"Total combinations with Decius Rules: {total_ic_sets}")


        
    

    '''
    Apply Decius Rules to possible Subgraphs
    '''    
    idx = 1
    combinations_sub_ics = []
    sub_ic_dicts = []
    fixed_bonds, fixed_angles, fixed_dihedrals = [], [], []
    n_r_fixed,n_phi_fixed,n_tau_fixed = [],[],[]

    print(ic_distribution.subgraph_ic_fixation(dis_subgraphs[0],submolecules[0]))
    for subgraph, submolecule in zip(dis_subgraphs,submolecules):
        bonds_sub, angles_sub, dihedrals_sub = ic_generator(subgraph)
        

        # TODO also fix this

        bonds_sub = [tuple(elem) for elem in bonds_sub]
        angles_sub = [tuple(elem) for elem in angles_sub]
        dihedral_sub = [tuple(elem) for elem in dihedrals_sub]


        fixed_bonds.extend(bonds_sub)
        writer.write_subheader(f"=== Spectral Bisection analysis subgraph {idx} ===")
        writer.write_line("Bonds")
        writer.write_data(bonds_sub)
        if len(angles_sub) != 0:
            writer.write_line("Angles")
            writer.write_data(angles_sub)
            fixed_angles.extend(angles_sub)
        if len(dihedrals_sub) != 0:
            writer.write_line("Dihedrals")
            writer.write_data(dihedrals_sub)
            fixed_dihedrals.extend(dihedrals_sub)

       
        # Apply Decius Rules to Subgraphs
        n_r_sub, n_phi_sub, n_tau_sub = general_nonlinear(submolecule,subgraph,bonds_sub)
        n_r_fixed.append(n_r_sub)
        n_phi_fixed.append(n_phi_sub)
        n_tau_fixed.append(n_tau_sub)

        writer.write_line("\n")
        writer.write_line("Application of Decius Rules to Submolecules")
        writer.write_line(f"Bond needed: {n_r_sub}")
        writer.write_line(f"Angles needed: {n_phi_sub}")
        writer.write_line(f"Dihedrals needed: {n_tau_sub}") 

        '''
        Build Submolecule IC sets as a Dictionary of Dictionareis
        '''



        
        bond_combs_sub = binomial_coefficient(len(bonds_sub),n_r_sub)
        if n_phi_sub > 0:
            angle_combs_sub = binomial_coefficient(len(angles_sub),n_phi_sub)
        else:
            angle_combs_sub = 1
        if n_tau_sub > 0: 
            dihedral_combs_sub = binomial_coefficient(len(dihedrals_sub),n_tau_sub)
        else:
            dihedral_combs_sub = 1

        writer.write_line(f"Bond subsets: {bond_combs_sub}")
        writer.write_line(f"Angle subsets: {angle_combs_sub}")
        writer.write_line(f"Dihedrals subsets: {dihedral_combs_sub}")

        total_ics_sub = bond_combs_sub * angle_combs_sub * dihedral_combs_sub 
        combinations_sub_ics.append(total_ics_sub)

                
        # Initialize the Nomodeco IC dict

        ic_dict_temp = {
            "bonds": [],
            "angles": [],
            "linear valence angles": [],
            "out of plane angles": [],
            "dihedrals": []
        
        }

        # initialize an empty dictionary with the total combination

        sub_ic_dict = {}


        ic_count = 0
        while ic_count < total_ics_sub:
            sub_ic_dict.setdefault(ic_count,ic_dict_temp.copy())
            ic_count +=1

        # Base case the combinations are 1
        if total_ics_sub == 1:
            for key in sub_ic_dict:
                sub_ic_dict[key]["bonds"].extend(bonds_sub)
                sub_ic_dict[key]["angles"].extend(angles_sub)
                sub_ic_dict[key]["dihedrals"].extend(dihedrals_sub)
        sub_ic_dicts.append(sub_ic_dict)
        
        # Then the case where the combinations are more then one
        if total_ics_sub > 1:
           # calculate the subsets as lists and then extend the lists
           bond_subsets = list(itertools.combinations(bonds_sub,n_r_sub))
           angle_subsets = list(itertools.combinations(angles_sub,n_phi_sub))
           dihedral_subsets = list(itertools.combinations(dihedrals_sub,n_tau_sub))
           
           # if the length of the subsets is 1 we can directly extend the lists to all the dictionaries
           if len(bond_subsets) == 1:
               for key in sub_ic_dict:
                   sub_ic_dict[key]["bonds"] = list(bond_subsets[0])
           if len(angle_subsets) == 1:
                for key in sub_ic_dict:
                   sub_ic_dict[key]["angles"] = list(angle_subsets[0]) 
           if len(dihedral_subsets) == 1:
                for key in sub_ic_dict:
                    sub_ic_dict[key]["dihedrals"] = list(dihedral_subsets[0])
           # if the length of the subset is more then one we append it modulo the total length
           if len(angle_subsets) > 1:
                for key in sub_ic_dict:
                     sub_ic_dict[key]["angles"] = list(angle_subsets[key % len(angle_subsets)])
           
           if len(dihedral_subsets) > 1:
               for key in sub_ic_dict:
                   sub_ic_dict[key]["dihedrals"] = list(dihedral_subsets[key % len(dihedral_subsets)]) # this zero here is a little bit confusing  

        writer.write_line(f"Total combination with Decius Rules: {total_ics_sub}")      
    

        idx += 1

    # Calulate how much coordinates are fixed in total
    n_r_fixed = sum(n_r_fixed)
    n_phi_fixed = sum(n_phi_fixed)
    n_tau_fixed = sum(n_tau_fixed)


    def combine_submolecule_ics(list_of_dicts):
        """
        Combines two submolecule ICs dicts taking into account all possible combinatoric cases
        """
        
        
        # basis case covered with the cartesian product
        if len(list_of_dicts[0]) == 1 and len(list_of_dicts[1]) == 1:
            combined = list_of_dicts[0]
            for key in combined[0]:
                combined[0][key].extend(list_of_dicts[1][0][key])

        # case one of the two dictionarys has length 1

        elif len(list_of_dicts[0]) == 1 and len(list_of_dicts[1]) > 1:
            # the total amount of combinations is given by the length of the second dictionary
            combinations = len(list_of_dicts[1])
            combined = list_of_dicts[1]
            # now loop over this combined dictionary and extend each one with the first dictionary
            for key in combined:
                combined[key]["bonds"] = list_of_dicts[0][0]["bonds"] + combined[key]["bonds"]  
                combined[key]["angles"] = list_of_dicts[0][0]["angles"] + combined[key]["angles"]
                combined[key]["dihedrals"] = list_of_dicts[0][0]["dihedrals"] + combined[key]["dihedrals"]
         
        elif len(list_of_dicts[1]) > 1 and len(list_of_dicts[0]) == 1:
            print("Hello")
        elif len(list_of_dicts[0]) > 1   and len(list_of_dicts[1]) > 1:
            # combinations are given by the multiplication of the lengths of the two dictionaries
            combinations = len(list_of_dicts[0]) * len(list_of_dicts[1])
            # make a copy of the first dictionary modulu the length of combinations
            combined = {}
            for i in range(combinations):
                combined.setdefault(i,list_of_dicts[0][i % len(list_of_dicts[0])].copy())
            # now append the second dictionary to the first modulu length of the total combinations
            for i in range(combinations):
                combined[i]["bonds"] = combined[i]["bonds"] + list_of_dicts[1][i % len(list_of_dicts[1])]["bonds"]
                combined[i]["angles"] = combined[i]["angles"] + list_of_dicts[1][i % len(list_of_dicts[1])]["angles"]
                combined[i]["dihedrals"] = combined[i]["dihedrals"] + list_of_dicts[1][i % len(list_of_dicts[1])]["dihedrals"]
        return combined
        
    
    combined_set = combine_submolecule_ics(sub_ic_dicts)
    


    writer.write_subheader("=== Combinations after Spectral Bisection ===")
    combinations_sub_ics = math.prod(combinations_sub_ics)
    writer.write_line(f"Possible combinations submolecules: {combinations_sub_ics}")
    # do this with sets because of the ordering issues
    fixed_bonds = [set(elem) for elem in fixed_bonds]
    fixed_angles = [set(elem) for elem in fixed_angles]
    fixed_dihedrals = [set(elem) for elem in fixed_dihedrals]
    # Calculate differences in sets, then Combinatorics
   
    bonds_choose = []
    for bond in bonds:
        if set(bond) not in fixed_bonds:
            bonds_choose.append(bond)
    angles_choose = []
    for angle in angles:
        if set(angle) not in fixed_angles:
            angles_choose.append(angle)
    dihedrals_choose = []
    for dihedral in dihedrals:
        if set(dihedral) not in fixed_dihedrals:
            dihedrals_choose.append(dihedral)
    

    # Now we write the function that builds all the possible combination

    def distribute_elements(ic_set, elements, n_choose, key):
        
        print(f"Combinatoric Distribution of Elements for {key}")
        print("=====================================")
        print("Length of IC Set: ", len(ic_set))
        print("Lenght of Elements: ", len(elements))
        print("Number of Elements to Choose: ", n_choose)

        # First we calculate the combinations
        combs = binomial_coefficient(len(elements),n_choose)
        # Now we calculate the total number of combinations
        total_combinations = len(ic_set) * combs

        # Make total_combinations number of copies of the ic_set
        result = {}
        for i in range(total_combinations):
            result.setdefault(i,ic_set[i % len(ic_set)].copy())
 
        
        # Now using itertools generate a list of all possible combinations, then loop

        comb_list = list(itertools.combinations(elements,n_choose))
        for i in range(total_combinations):
            result[i][key] = result[i][key] + list(comb_list[i % len(comb_list)])
        
        # Print the result
        print("=====================================")
        print(f"New IC_Dict has lenght: {len(result)}")

        return result

    ic_dict = distribute_elements(combined_set,bonds_choose, (n_r-n_r_fixed),"bonds")
    ic_dict = distribute_elements(ic_dict,angles_choose, (n_phi-n_phi_fixed),"angles")
    ic_dict = distribute_elements(ic_dict,dihedrals_choose, (n_tau-n_tau_fixed),"dihedrals")


    bonds_choose_combs = binomial_coefficient(len(bonds_choose), (n_r - n_r_fixed))
    angles_choose_combs = binomial_coefficient(len(angles_choose), (n_phi - n_phi_fixed))
    dihedrals_choose_combs = binomial_coefficient(len(dihedrals_choose),(n_tau - n_tau_fixed))

    writer.write_line(f"Bond subsets: {bonds_choose_combs}")
    writer.write_line(f"Angle subsets: {angles_choose_combs}")
    writer.write_line(f"Dihedral subsets: {dihedrals_choose_combs}")

    spectral_bisection_sets = combinations_sub_ics * bonds_choose_combs * angles_choose_combs * dihedrals_choose_combs
    writer.write_line(f"Spectral Bisection Total Sets: {spectral_bisection_sets}")

    # Write a seperate outpufile with the ic_dict
    with open("ic_dict_output.txt", "w") as outfile:
        for key, value in ic_dict.items():
            outfile.write(f"{key}: {value}\n")

    
    '''
    Graph Pruning approach

    Take the pruned graph, fix all ICs for the pruned graph and append missing coordinates
    '''
    writer.write_subheader("=== Graph Pruning Approach ===")

    # Remove this maybe later
    print("=====================================")
    print("Graph Pruning Approach")
    print("=====================================")

    # First of we need to generate the ICs for the pruned graph
    pruned_bonds, pruned_angles, pruned_dihedrals = ic_generator(pruned_subgraph)
    # calculate the Number of ICs needed (DECIUS)
    n_r_pruned, n_phi_pruned, n_tau_pruned = general_nonlinear(pruned_submolecules,pruned_subgraph,pruned_bonds) 

    # Initializte the IC dict for the pruned graph
    ic_dict_pruned = {
        "bonds": [],
        "angles": [],
        "linear valence angles": [],
        "out of plane angles": [],
        "dihedrals": []
    }
    # Make one initial entry in the big ic dict and then use the distribute elements function
    ic_dict_pruned = {0: ic_dict_pruned}
    ic_dict_pruned = distribute_elements(ic_dict_pruned,pruned_bonds, n_r_pruned,"bonds")
    ic_dict_pruned = distribute_elements(ic_dict_pruned,pruned_angles,n_phi_pruned,"angles")
    if n_tau_pruned > 0:
        ic_dict_pruned = distribute_elements(ic_dict_pruned,pruned_dihedrals,n_tau_pruned,"dihedrals")
    
    # evaluate the fixed bonds
    fixed_bonds_pruned = [set(elem) for elem in pruned_bonds]
    fixed_angles_pruned = [set(elem) for elem in pruned_angles]
    fixed_dihedrals_pruned = [set(elem) for elem in pruned_dihedrals]
    
    # now calculate the differences
    bonds_choose_pruned = []
    for bond in bonds:
        if set(bond) not in fixed_bonds_pruned:
            bonds_choose_pruned.append(bond)
    angles_choose_pruned = []
    for angle in angles:
        if set(angle) not in fixed_angles_pruned:
            angles_choose_pruned.append(angle)
    dihedrals_choose_pruned = []
    for dihedral in dihedrals:
        if set(dihedral) not in fixed_dihedrals_pruned:
            dihedrals_choose_pruned.append(dihedral)
    
    # Now again use our distribute elements function

    ic_dict_pruned = distribute_elements(ic_dict_pruned,bonds_choose_pruned, (n_r - n_r_pruned),"bonds")
    ic_dict_pruned = distribute_elements(ic_dict_pruned,angles_choose_pruned, (n_phi - n_phi_pruned),"angles")
    if n_tau - n_tau_pruned > 0:
        ic_dict_pruned = distribute_elements(ic_dict_pruned,dihedrals_choose_pruned, (n_tau - n_tau_pruned),"dihedrals")


    # Write decius numbers to out

    writer.write_line(f"Bonds needed: {n_r_pruned}")
    writer.write_line(f"Angles needed: {n_phi_pruned}")
    writer.write_line(f"Dihedrals needed: {n_tau_pruned}")

    # Write the total number of combinations
    writer.write_line(f"Total IC sets: {len(ic_dict_pruned)}")

    # Write a seperate outpufile with the ic_dict
    with open("ic_dict_output_pruned.txt", "w") as outfile:
        for key, value in ic_dict_pruned.items():
            outfile.write(f"{key}: {value}\n")

    '''
    Kerighan_Lin Algorithm
    ''' 

    kerighan_lin = graph_analysis.kerighan_lin_algorithm(G)
     
    # make into list of lists
    kerighan_lin = [list(elem) for elem in kerighan_lin]
    print(kerighan_lin)




if __name__ == "__main__":
    main()