####################################################################################################
# This module contains functions to visualize the graph.
####################################################################################################


import matplotlib.pyplot as plt
import networkx as nx
import numpy as np 
from matplotlib.lines import Line2D


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

    # Create custom legend
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'Degree {deg}', 
                              markerfacecolor=color, markersize=10) 
                       for deg, color in degree_colors.items()]
    fig.legend(handles=legend_elements, loc='upper right', title="Node Degree")


    plt.tight_layout()
    plt.savefig(f"{figname}.png", dpi=600)
    
    #Do show if you want nice animation
    #plt.show()