from itertools import product

from pyqubo import Qbit


def graph_to_qbits(G):
    """
    Get a hash table of string variable labels -> Qbit variables.

    :param G: Networkx.Graph
    :return: copy of G, dict of {node label: Qbit}
    """
    # copy so we don't affect the original graph
    G_copy = G.copy()

    # Create variable label for nodes
    for u in sorted(G_copy.nodes):
        G_copy.node[u]["label"] = "x_{}".format(u)

    return G_copy, {uvar['label']: Qbit(uvar['label']) for _, uvar in G_copy.nodes(data=True)}


def graph_to_qbits_twovar(G):
    """
    For coloring problems that require a choice of 'color' for each node. Thus, we get 'x_i_j',
        where i is the node id, and j is the 'color' choice for node i.

    :param G: NetworkX.Graph
    :return: G_copy,  dict of {: Qbit}
    """
    # Create variables from nodes and colors
    # x_v_i = 1 iff vertex v has color i, zero otherwise.
    node_color_str = ["x_{}_{}".format(u[0], u[1]) for u in product(G.nodes, range(num_colors))]
