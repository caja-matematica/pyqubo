# -*- coding: utf-8 -*-
from sympy import flatten
from itertools import combinations
from math import log

from pyqubo.utils import graph_converter as gc

from qubo.core.symbolic_utils import *
from qubo.core.qubo_response import QuboResponse


def generate_vertex_cover_qubo(G, A=10, B=5):
    """
    Lucas, 4.3 Vertex Cover

    Given an undirected graph $G= (V, E)$, what is the smallest number of vertices
    that can be "colored" such that every edge is incident to a colored vertex.

    Choose B < A, as if we uncolor any vertex that ruins the solution, at least one edge will
    no longer connect to a colored vertex.

    Args:
        G: NetworkX Graph
        A: float
            positive scalar factor
        B: float
            positive scalar factor (see documentation above for ratio of A to B)

    Returns:
         :class:`.qubo_response.QuboResponse`
    """

    if not H:
        raise ValueError("Cannot handle empty graph")

    # G is a copy of H in which node attributes have been added
    G_copy, node_hash = gc.graph_to_qbits(G)

    edge_constraint = []
    for u,v in G_copy.edges:
        x_u = node_hash[G_copy.node[u]['label']]
        x_v = node_hash[G_copy.node[v]['label']]
        edge_constraint.append((1 - x_u) * (1 - x_v))

    HA = A * sum(edge_constraint)
    HB = B * sum(node_hash.values())
    H = HA + HB

    return H
    # return QuboResponse(H, symbol_hash)


def generate_vertex_cover_qubo2(G, A=10, B=5):
    """
    Lucas, 4.3 Vertex Cover

    Given an undirected graph $G= (V, E)$, what is the smallest number of vertices
    that can be "colored" such that every edge is incident to a colored vertex.

    Choose B < A, as if we uncolor any vertex that ruins the solution, at least one edge will
    no longer connect to a colored vertex.

    Args:
        G: NetworkX Graph
        A: float
            positive scalar factor
        B: float
            positive scalar factor (see documentation above for ratio of A to B)

    Returns:
         :class:`.qubo_response.QuboResponse`
    """

    if not G:
        raise ValueError("Cannot handle empty graph")

    # Create variables from nodes and get edgelist
    node_str = ["x_{}".format(u) for u in sorted(G.nodes)]

    # QSymbol defined below (for now)
    node_vars = [QSymbol(x) for x in node_str]
    node_hash = {v.get_subscript(): v for v in node_vars}

    # enumeration of variables to place in a matrix
    symbol_hash = {sym: idx for idx, sym in enumerate(flatten(node_vars))}

    edge_constraint = []
    for e in G.edges:
        x_u = node_hash[e[0]]
        x_v = node_hash[e[1]]
        edge_constraint.append((1 - x_u) * (1 - x_v))

    # sum over all edges to get that each edge has at least one colored vertex
    HA = A * sum(edge_constraint)

    # minimize number of colored vertices
    HB = B * sum(node_hash.values())

    H = HA + HB

    return QuboResponse(H, symbol_hash)


def generate_graph_color_qubo(G, num_colors=4, A=1):
    """
    6.1 Graph Coloring

    Given an undirected graph $G = (V, E)$, and a set of n colors, is it possible to color each vertex in the
    graph with a specific color, such that no edge connects two vertices of the same color?


    Args:
        G: NetworkX Graph
        A: float
            positive scalar factor
        num_colors: int
            number of cliques to form

    Returns:
         :class:`.qubo_response.QuboResponse`
    """

    if not G:
        raise ValueError("Cannot handle empty graph")

    # Create variables from nodes and colors
    # x_v_i = 1 iff vertex v has color i, zero otherwise.
    node_color_str = ["x_{}_{}".format(u[0], u[1]) for u in product(G.nodes, range(num_colors))]

    G_copy, node_hash = gc.graph_to_qbits(H)

    # convert string vars to Sympy Symbols
    # node_color_vars = [TwoSubscriptQSymbol(x) for x in node_color_str]

    # puts i,j subscripts into easily searchable hash table
    # eg., color_options_hash = {(i,j): symbol}
    # and symbol_hash = {symbol: (i,j)}
    # color_options_hash = {v.get_subscript(): v for v in node_color_vars}
    symbol_hash = {sym: idx for idx, sym in enumerate(node_color_vars)}

    node_groups = groupby_index(node_hash, idx=0)

    node_constraints = []
    for _, group in node_groups:
        # fix i, result = \sum_j {x_i_j}
        node_const = sum_over_fixed_first_subscript(color_options_hash, group)
        # (1 - \sum_j {x_i_j} )^2
        node_constraints.append((1 - node_const) ** 2)

    # sum over each vertex constraint
    H1 = sum(node_constraints)

    edge_penalties = []
    for e in G.edges:
        x_u_all = get_edge_colors(e, node_hash, node_id=0)
        x_v_all = get_edge_colors(e, color_options_hash, node_id=1)
        # sum_i {x_u,i * x_v,i}
        edge_uv_sum = sum(map(lambda z: z[0] * z[1], zip(x_u_all, x_v_all)))
        edge_penalties.append(edge_uv_sum)

    # sum over all the edges
    H2 = sum(edge_penalties)

    H = A * (H1 + H2)

    return QuboResponse(H, symbol_hash)


def generate_graph_color_qubo2(G, num_colors=4, A=1):
    """
    6.1 Graph Coloring

    Given an undirected graph $G = (V, E)$, and a set of n colors, is it possible to color each vertex in the
    graph with a specific color, such that no edge connects two vertices of the same color?


    Args:
        G: NetworkX Graph
        A: float
            positive scalar factor
        num_colors: int
            number of cliques to form

    Returns:
         :class:`.qubo_response.QuboResponse`
    """

    if not G:
        raise ValueError("Cannot handle empty graph")

    # Create variables from nodes and colors
    # x_v_i = 1 iff vertex v has color i, zero otherwise.
    node_color_str = ["x_{}_{}".format(u[0], u[1]) for u in product(G.nodes, range(num_colors))]

    # convert string vars to Sympy Symbols
    node_color_vars = [TwoSubscriptQSymbol(x) for x in node_color_str]

    # puts i,j subscripts into easily searchable hash table
    # eg., color_options_hash = {(i,j): symbol}
    # and symbol_hash = {symbol: (i,j)}
    color_options_hash = {v.get_subscript(): v for v in node_color_vars}
    symbol_hash = {sym: idx for idx, sym in enumerate(node_color_vars)}

    node_groups = groupby_index(color_options_hash, idx=0)

    node_constraints = []
    for _, group in node_groups:
        # fix i, result = \sum_j {x_i_j}
        node_const = sum_over_fixed_first_subscript(color_options_hash, group)
        # (1 - \sum_j {x_i_j} )^2
        node_constraints.append((1 - node_const) ** 2)

    # sum over each vertex constraint
    H1 = sum(node_constraints)

    edge_penalties = []
    for e in G.edges:
        x_u_all = get_edge_colors(e, color_options_hash, node_id=0)
        x_v_all = get_edge_colors(e, color_options_hash, node_id=1)
        # sum_i {x_u,i * x_v,i}
        edge_uv_sum = sum(map(lambda z: z[0] * z[1], zip(x_u_all, x_v_all)))
        edge_penalties.append(edge_uv_sum)

    # sum over all the edges
    H2 = sum(edge_penalties)

    H = A * (H1 + H2)

    return QuboResponse(H, symbol_hash)


def generate_clique_cover_qubo(G, A=1, B=1, num_cliques=4):

    """
    6.2 Clique Cover

    The clique cover problem, for an undirected graph $G = (V, E)$, is the following: given $n$ colors, we assign
    a distinct color to each vertex of the graph. Let $W_1,\ldots, W_n$ be the subsets of $V$ corresponding to each
    color, and $E_{W_1},\ldots, E_{W_n}$ the edge set restricted to edges between vertices in the Wi sets.
    The clique cover problem asks whether or not $(W_i, E_{W_i})$ is a complete graph for each $W_i$
    (i.e., does each set of colored vertices form a clique?).


    There is a new and important function here: `get_two_subscript_symbols`.
    Basically, this problem requires that for eac vertex $v$ in $G$, we have to account for `num_variations`
    (`num_colors` in this case) options for each $v$. Hence, instead of creating the
    variable $x_{v}$ for the node $v$, if `num_colors` = 4 we create $x_{v,0}, x_{v,1}, x_{v,2}, x_{v,3}$.
    `get_two_subscript_symbols` returns the corresponding $(u,j) \rightarrow x_{u,j}$ mapping,
    plus the symbol to index mapping for the QUBO matrix.

    See Lucas, Section 2.3 for a detailed discussion of the ratio of A/B. In short, choosing

        A = (\Delta +2)B,

    where \Delta = maximal degree of G.

    Args:
        G: NetworkX Graph
        A: float
            positive scalar factor
        B: float
            positive scalar factor (see documentation above for ratio of A to B)
        num_cliques: int
            number of cliques to form

    Returns:
         :class:`.qubo_response.QuboResponse`
    """

    if not G:
        raise ValueError("Cannot handle empty graph")

    node_hash, symbol_hash = get_two_subscript_symbols(G.nodes, num_variations=num_cliques)

    # This will group variables by the first index, eg., v in x_{v,i}
    node_groups = groupby_index(node_hash)

    node_constraints = []
    for _, group in node_groups:
        node_const = sum_over_fixed_first_subscript(node_hash, group)
        node_constraints.append((1 - node_const) ** 2)

    # sum over each vertex constraint
    HA = A * sum(node_constraints)

    # This loop fies the second subscript (the "variations" subscript, or "color" in this case).
    # Fix color terms (second indices) and constrain based on allowable colors per node
    color_sum_terms = []
    for color in range(num_cliques):
        color_constraint_sum = sum_over_fixed_second_subscript(node_hash, second_idx=color)
        HB1 = 0.5 * (color_constraint_sum - 1) * color_constraint_sum

        # This is a sum over all edges (u,v), with color fixed as 'c'. We get \sum_{u,v} x_{u,c}*x_{v,c}
        HB2 = sum_edges_over_fixed_second_subscript(G.edges, node_hash, color)

        color_sum_terms.append(HB1 - HB2)

    H = HA + B * sum(color_sum_terms)

    return QuboResponse(H, symbol_hash)


def generate_hamiltonian_cycle_qubo(G, A=1, exclude_diag=False):
    """
    Let $G = (V, E)$, and $N = |V |$. The graph can either be directed or undirected; our method of solution
    will not change. The Hamiltonian path problem is as follows: starting at some node in the graph, can
    one travel along an edge, visiting other nodes in the graph, such that one can reach every single node in
    the graph without ever returning to the same node twice? The Hamiltonian cycles problem asks that, in
    addition, the traveler can return to the starting point from the last node he visits.

    Args:
        G: NetworkX Graph
        A: float
            positive scalar factor
        exclude_diag: bool
            Controls whether to allow diagonal elements to explicitly exclude self-loops for us in
            Traveling Salesman problems. For instance, the 'two-visit' formulation of the traveling salesman problem
            leverages explicit exclusion of the diagonal by setting up constraints such as

               (ab + ac + ad + ae + af + ag - 2)^2,

            where a, b, ..., g are cities, and each one much be visited two, once coming and once going. Note, though, that
            aa was excluded. Alternatively, we can enforce such behavior directly through the structure of the energy
            function. See Lucas, 7.2 and 7.3.

    Returns:
         :class:`.qubo_response.QuboResponse`
    """
    if not G:
        raise ValueError("Cannot handle empty graph")

    node_options_hash, symbol_hash = get_two_subscript_symbols(G.nodes, num_variations=len(G.nodes))

    # group nodes my vertex
    node_groups = groupby_index(node_options_hash, idx=0)

    # every vertex must appear in a cycle
    node_constraints = []
    for _, group in node_groups:
        node_const = sum_over_fixed_first_subscript(node_options_hash, group)
        node_constraints.append((1 - node_const) ** 2)

    H1 = sum(node_constraints)

    # there must be a jth node in each cycle
    label_constraints = []
    for j in range(len(G.nodes)):
        label_const = sum_over_fixed_second_subscript(node_options_hash, second_idx=j)
        label_constraints.append((1 - label_const) ** 2)

    H2 = sum(label_constraints)

    # energy penalty if path uv chosen but does not exist in G
    possible_edges = set(combinations(G.nodes, 2))
    non_edges = possible_edges.difference(G.edges)
    num_labels = len(G.nodes)

    # edge penalty defined locally
    H3 = sum([edge_penalty(e, num_labels, node_options_hash) for e in non_edges])

    H = A * (H1 + H2 + H3).expand()

    return QuboResponse(H, symbol_hash)


def generate_tsp_qubo(G, path_distances, A=1, B=1):
    """
    7.2 Traveling Salesman (Lucas version)

    The traveling salesman problem for a graph $G = (V, E)$, where each edge $uv$ in the graph has a weight
    $W_{uv}$ associated to it, is to find the Hamiltonian cycle such that the sum of the weights of each edge in the
    cycle is minimized. Typically, the traveling salesman problem assumes a complete graph, but we have the
    technology developed to solve it on a more arbitrary graph.

    Given a file with cities and distance, see converters.distance_file_graph_converter() for sample
    methodology for converting a distance file to a Graph and a path_distance dict.

    Possible heuristic for A and B: 0 < B * max(W_{uv}) < A

    Args:
        G: NetworkX Graph
        path_distances (dict): {(city_a, city_b): distance}.
            Eg., {(0, 1): 2230, (0, 2): 1631, (0, 3): 1566, ... (4, 6): 916, (5, 6): 702}
        A (positive float): scalar factor
        B (positive float): scalar factor (see above for relationship of A to B)

    Returns:
         :class:`.qubo_response.QuboResponse`
    """
    # we exclude the default symbol hash since in this case we want "stay off of the diagonal"
    node_options_hash, symbol_hash = get_two_subscript_symbols(G.nodes, num_variations=len(G.nodes))

    quboHA = generate_hamiltonian_cycle_qubo(G, A=A)
    HA = quboHA.get_hamiltonian()

    # distance minimization
    HB = B * sum([path_distances[e] * edge_penalty(e, len(G.nodes), node_options_hash) for e in G.edges])

    H = HA + HB

    return QuboResponse(H, symbol_hash)


def generate_tsp_qubo_legacy(distance_file, city_hash, num_cities=7, gamma=1):
    """

    Args:
        distance_file: str
            path to file containing distances between cities (typically a ".b" file)
        city_hash: dict
            string name -> integer hash. Eg., {"a": 0, "b": 1, ..., "g": 6}
        num_cities: int
            Number of cities (TODO: we can get this from the hash)
        gamma: float
            Multiplier to align scale of constraints to distance scale

    Returns:
         :class:`.qubo_response.QuboResponse`
    """

    with open(distance_file) as fh:
        dists = fh.readlines()

    # Build up possible paths == variables. A path a -> b is represented by a variable x_{a}_{b}
    # TODO: we no longer enforce upper triangular, so we need to create a new function to handle that necessity here
    node_hash, symbol_hash = get_two_subscript_symbols(range(num_cities), num_variations=num_cities)

    path_dists = {}
    for line in dists:
        line = line.strip()
        if not line:
            continue
        if not line.startswith("D"):
            continue
        two_cities, d = line.split(" = ")

        # strip off the "D" and split the path into start/end
        start, end = tuple(two_cities[1:])
        start_num = city_hash[start]
        end_num = city_hash[end]

        path_dists[(start_num, end_num)] = int(d)

    dterm = []
    for city_path, city_var in node_hash.items():
        try:
            dterm.append(path_dists[city_path] * city_var)
        # we only keep i < j due to symmetric distances
        # Key error thrown when path_{j,i} not found in paths dict
        except KeyError:
            pass

    H1 = sum(dterm)

    first_city_grps = groupby_index(node_hash, idx=0)

    starting_city_constraints = []

    # \sum {x_a_b} - 2 => each city is visited twice, once entering, once exiting, in terms of path variables.
    for _, path_gp in first_city_grps:
        starting_city_constraints.append((sum([node_hash[path] for path in path_gp]) - 2) ** 2)

    H2 = gamma * sum(starting_city_constraints)

    H = H1 + H2

    return QuboResponse(H, symbol_hash)


def generate_number_partitioning_qubo(list_to_partition):

    """
    Lucas, 2.1 Number Partitioning

    Given a set of N positive numbers $S = {n_1, n_2... n_N}$,
    is there a partition of this set of numbers into two disjoint subsets
    R and S - R, so that the sum of the elements in both subsets is the
    same

    Lucas's formulation is for the Ising model. This formulation is for
    QUBO, obtained by standard variable substitution.

    Args:
        list_to_partition: List

    Returns:
         :class:`.qubo_response.QuboResponse`

    """

    if len(list_to_partition) == 0:
        raise ValueError("Cannot handle empty list")

    # Create variables from nodes and get edgelist
    node_str = ["x_{}".format(u) for u in range(len(list_to_partition))]

    # QSymbol defined below (for now)
    node_vars = [QSymbol(x) for x in node_str]
    node_hash = {v.get_subscript(): v for v in node_vars}

    # enumeration of variables to place in a matrix
    symbol_hash = {sym: idx for idx, sym in enumerate(flatten(node_vars))}

    sum_all_numbers = sum(list_to_partition)

    # summation of one partition
    sum_temp = [2 * list_to_partition[x] * node_hash[x] for x in range(len(list_to_partition))]
    sum_partition = sum(sum_temp)

    # The QUBO is the difference between the total sum (both partitions) and
    # one of the partition. We want to minimize this difference.
    H = (sum_all_numbers - sum_partition) * (sum_all_numbers - sum_partition)

    return QuboResponse(H, symbol_hash)


def generate_maxcut_qubo(G):

    """
    The Maximum Cut problem

    Mark Lewis and Fred Glover, https://arxiv.org/pdf/1705.09844.pdf

    Given a graph $G= (V, E)$, one wants a subset S of the vertex set such
    that the number of edges between S and the complementary subset is as
    large as possible

    Args:
        G: NetworkX Graph

    Returns:
         :class:`.qubo_response.QuboResponse`

    """
    if not G:
        raise ValueError("Cannot handle empty graph")

    # Create variables from nodes and get edgelist
    node_str = ["x_{}".format(u) for u in sorted(G.nodes)]

    # QSymbol defined below (for now)
    node_vars = [QSymbol(x) for x in node_str]
    node_hash = {v.get_subscript(): v for v in node_vars}

    # enumeration of variables to place in a matrix
    symbol_hash = {sym: idx for idx, sym in enumerate(flatten(node_vars))}

    edge_constraint = []
    for e in G.edges:
        x_u = node_hash[e[0]]
        x_v = node_hash[e[1]]
        edge_constraint.append(2 * x_u * x_v - x_u - x_v)

    H = sum(edge_constraint)

    return QuboResponse(H, symbol_hash)


def generate_knapsack_qubo(valuez, weightz, W, A=10, B=1):

    """
    Lucas, 5.2 Knapsack with Integer Weights

    We have a list of N objects (weightz), labeled by indices i, with
    the weight of each object given by w_i, and each has a value given by
    c_i. The knapsack can only hold weight W. The goal is to maximize the
    total value C = sum_i c_i, while filling the backpack as close to W
    as possible.

    A controls the constraints; B controls the costs. The goal is to find
    a suitable tradeoff between constraints and costs so that the overall
    goal is met.

    Note: this algorithm, in Lucas's paper, does not converge satisfactorily.
    Lucas seems to have updated it. Please see "knapsack2," below.

    Args:
        valuez: List of costs
        weightz: List of weightz
        W: maximum weight allowed in backpack
        A: float
            positive scalar factor
        B: float
            positive scalar factor

    Returns:
         :class:`.qubo_response.QuboResponse`

    """

    if not valuez or not weightz:
        raise ValueError("Cannot handle empty lists.")

    # Creating string variables for x (determining whether the object is in
    # the knapsack), and y (whether the final weight is n)
    node_x_str = ["x_{}".format(u) for u in range(len(weightz))]
    node_y_str = ["y_{}".format(u) for u in range(W)]

    # Creating symbol variables for x and y
    node_x_vars = [QSymbol(x) for x in node_x_str]
    node_y_vars = [QSymbol(y) for y in node_y_str]

    # Create one array of all variables - first x and then y
    node_vars = node_x_vars + node_y_vars

    # Form lookup table for the indices of the individual variables
    symbol_hash = {sym: idx for idx, sym in enumerate(flatten(node_vars))}

    # Set up the sum of the y variables
    sum_y = sum(node_y_vars)

    # Set up the sum of n * the y variables, except note that it has to
    # start from 1.
    ny_val = [n * node_y_vars[n] for n in range(len(node_y_vars))]

    # Set up the wx term
    wx_val = [weightz[k] * node_x_vars[k] for k in range(len(weightz))]

    # Form Lucas's first sum, Equation (49)
    HA = ((1 - sum_y) ** 2) + ((sum(ny_val) - sum(wx_val)) ** 2)

    # Form the cx sum, for Lucas's second sum, Equation (50)
    cx_val = [valuez[k] * node_x_vars[k] for k in range(len(valuez))]
    HB = -B * sum(cx_val)

    H = (A * HA) + HB

    return QuboResponse(H, symbol_hash)


def generate_knapsack2_qubo(valuez, weightz, W, A=10):

    """
    Lucas, 5.2 Knapsack with Integer Weights - Second Algorithm

    We have a list of N objects (weightz), labeled by indices i, with
    the weight of each object given by w_i, and each has a value given by
    c_i. The knapsack can only hold weight W. The goal is to maximize the
    total value C = sum_i c_i, while filling the backpack as close to W
    as possible.

    A controls the constraints; B controls the costs. The goal is to find
    a suitable tradeoff between constraints and costs so that the overall
    goal is met.

    This algorithm was found on Lucas's slides from 2014. It converges
    much more readily than "knapsack".

    Args:
        valuez: List of costs
        weightz: List of weightz
        W: maximum weight allowed in backpack
        A: float
            positive scalar factor

    Returns:
         :class:`.qubo_response.QuboResponse`

    """

    if not valuez or not weightz:
        raise ValueError("Cannot handle empty lists.")

    logw = int(log(W)) + 1
    m = logw - 1

    node_x_str = ["x_{}".format(u) for u in range(len(weightz))]
    node_y_str = ["y_{}".format(u) for u in range(logw)]

    node_x_vars = [QSymbol(x) for x in node_x_str]
    node_y_vars = [QSymbol(y) for y in node_y_str]

    node_vars = node_x_vars + node_y_vars
    symbol_hash = {sym: idx for idx, sym in enumerate(flatten(node_vars))}

    sumy = (W + 1 - (2 ** m)) * node_y_vars[m]
    for n in range(m):
        sumy += (2 ** n) * node_y_vars[n]

    # Set up the wx term
    wx_val = [weightz[k] * node_x_vars[k] for k in range(len(weightz))]

    HA = (sumy - sum(wx_val)) ** 2

    cx_val = [valuez[k] * node_x_vars[k] for k in range(len(valuez))]
    HB = -sum(cx_val)
    H = (A * HA) + HB

    return QuboResponse(H, symbol_hash)


def generate_maximum_independent_set_qubo(G, A=10):

    """
    The Maximum Independent Set problem

    Mark Lewis and Fred Glover, https://arxiv.org/pdf/1705.09844.pdf

    Args:
        G: NetworkX Graph
        A: float
            positive scalar factor

    Returns:
         :class:`.qubo_response.QuboResponse`

    """
    if not G:
        raise ValueError("Cannot handle empty graph")

    # Create variables from nodes and get edgelist
    node_str = ["x_{}".format(u) for u in sorted(G.nodes)]

    # QSymbol defined below (for now)
    node_vars = [QSymbol(x) for x in node_str]
    node_hash = {v.get_subscript(): v for v in node_vars}

    # enumeration of variables to place in a matrix
    symbol_hash = {sym: idx for idx, sym in enumerate(flatten(node_vars))}

    edge_constraint = []
    for e in G.edges:
        x_u = node_hash[e[0]]
        x_v = node_hash[e[1]]
        edge_constraint.append(x_u * x_v)

    # sum over all edges to get that each edge has at least one colored vertex
    HA = A * sum(edge_constraint)

    # minimize number of colored vertices
    HB = -sum(node_hash.values())

    H = HA + HB

    return QuboResponse(H, symbol_hash)


def generate_job_assignment_qubo(weightz, A=1, B=1):

    """
    The Job Assignment Problem

    For example,
    http://www.math.harvard.edu/archive/20_spring_05/handouts/assignment_overheads.pdf

    Given a square matrix $W(i,j)$, containing how long each worker i needs to
    finish job j, find the optimal assignment of the jobs to workers, to
    minimize the total time invested.

    Args:
        Weightz: Square matrix of job weights
        A: float
            positive scalar factor
        B: float
            positive scalar factor

    Returns:
         :class:`.qubo_response.QuboResponse`

    """
    if not weightz:
        raise ValueError("Cannot handle empty weight matrix")

    lenw = len(weightz)
    node_str = ["x_{}_{}".format(u[0], u[1]) for u in product(range(lenw), range(lenw))]
    node_vars = [TwoSubscriptQSymbol(x) for x in node_str]
    node_hash = {v.get_subscript(): v for v in node_vars}
    symbol_hash = {sym: idx for idx, sym in enumerate(flatten(node_vars))}

    cost_matrix = 0
    for i in range(lenw):
        for j in range(lenw):
            cost_matrix += weightz[i][j] * node_hash[(i, j)]

    HA = 0
    for i in range(lenw):
        row_constraint = 1
        for j in range(lenw):
            row_constraint -= node_hash[(i, j)]
        HA += row_constraint ** 2
        col_constraint = 1
        for j in range(lenw):
            col_constraint -= node_hash[(j, i)]
        HA += col_constraint ** 2

    HB = B * cost_matrix

    H = (A * HA) + HB

    return QuboResponse(H, symbol_hash)
