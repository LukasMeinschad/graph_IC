####################################################################################################
# Description: This module contains general functions that are used by various other modules.
####################################################################################################

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from string import digits
from collections import Counter


# Math Functions

def binomial_coefficient(n,k):
    """ 
    Calculates the Binomial Coefficient
    """

    return math.comb(n,k)



def remove_digits(s):
    """
    Removes digits from a string.

    Args:
        s (str): The input string.

    Returns:
        str: The string with digits removed.
    """
    return s.translate(str.maketrans('', '', '0123456789'))




def chemical_formula(molecule):
    """
    Extracts the Chemical Formula from a molecule
    """
    # First store the atoms in a list
    atoms = []
    for atom in molecule:
        atoms.append(atom.symbol)
    
    # remove the numbering
    atoms = [remove_digits(atom) for atom in atoms]

    # Count the elements
    element_counts = Counter(atoms)

    formula = ''.join(f"{element}{count if count > 1 else ''}" for element, count in sorted(element_counts.items()))

    return formula



     