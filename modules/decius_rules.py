############################################
# Decius Rules
############################################

import networkx as nx
import numpy as np



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
