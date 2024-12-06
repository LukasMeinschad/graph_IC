from mendeleev.fetch import fetch_table
import string
import numpy as np
import itertools
import pandas as pd

def get_bond_information():
    df = fetch_table("elements")
    bond_info = df.loc[:, ["symbol","covalent_radius_pyykko","vdw_radius"]]
    
    bond_info.set_index("symbol",inplace=True)
    bond_info /= 100
    return bond_info

bond_info = get_bond_information()

# Calculate the theoretical length

def theoretical_length(symbol1,symbol2):
    

    rad_a = bond_info.loc[symbol1.strip(string.digits)]
    rad_b = bond_info.loc[symbol2.strip(string.digits)]

    return rad_a[0] + rad_b[0]

def theoretical_length_vdw(symbol1,symbol2):

    rad_a = bond_info.loc[symbol1.strip(string.digits)]
    rad_b = bond_info.loc[symbol2.strip(string.digits)]

    return rad_a[1] + rad_b[1]

def actual_length(atom1,atom2):
    return abs(np.linalg.norm(np.array(atom1.coordinates) - np.array(atom2.coordinates)))

# Caculate the degree of covalence table 

def degree_of_covalance(molecule):
        """
        Calculates the degree of covalance between all the combinations of two atoms
        
        Reference  https://doi.org/10.1002/qua.21049

        Returns:
            a python dictionary with the two atoms as a tuple and the value of the degree of covalance
        """
        degofc = {} # initialize array
        for atom_a,atom_b in itertools.combinations(molecule,2):
            value = np.exp(-(abs(actual_length(atom_a,atom_b))/theoretical_length(atom_a.symbol,atom_b.symbol) - 1))
            degofc.update({(atom_a.symbol,atom_b.symbol) : value})
        return degofc



# Establish all covalent bonds
def covalent_bonds(molecule, degofc_table) -> list:
    
    combinations = [(atom_a.symbol,atom_b.symbol) for atom_a,atom_b in itertools.combinations(molecule,2)]
    hits = []
    for key,value in degofc_table.items():
        if key in combinations and value > 0.75:
            hits.append(key)
    return hits