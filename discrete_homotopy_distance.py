import copy
from sympy import Poly, symbols
from sympy.abc import x
from copy import deepcopy
import math
import networkx as nx

import random


import matplotlib.pyplot as plt
import time
from matplotlib.ticker import MaxNLocator

import plotly.graph_objects as go

from collections import deque, defaultdict


def succs(h_poly):
    """
    h_poly: initial homotopy polynomial (dict of {exponent: coefficient})
    or any input convertible to Poly within the sympy library e.g. sympy.Poly

    Returns:
        (successors, succs_list)
        successors: list of sympy.Poly polynomials
        succs_list: list of dictionaries ({exponent: coefficient})
    """

    h_poly = Poly(h_poly, x)
    monoms = h_poly.monoms()
    coeffs = h_poly.coeffs()
    alpha = dict(zip([m[0] for m in monoms], coeffs))

    successors = []
    succs_list = []
    seen = set()

    for y in list(alpha.keys()):
        g = copy.deepcopy(alpha)
        #g = dict(alpha)
        if g[y] > 0:
            g[y] -= 1

            for z in list(g.keys()):
                f = copy.deepcopy(g)
                #f = dict(g)
                if f[z] > 0:
                    f[z] -= 1
                    f[y + z] = f.get(y + z, 0) + 1

                    f_clean = {k: v for k, v in f.items() if v != 0}
                    f_tuple = tuple(sorted(f_clean.items()))
                    if f_tuple not in seen:
                        seen.add(f_tuple)
                        poly_f = sum(c * x**e for e, c in f_clean.items())
                        successors.append(Poly(poly_f, x))
                        succs_list.append(f_clean)

            g[y + 1] = g.get(y + 1, 0) + 1
            g_clean = {k: v for k, v in g.items() if v != 0}
            g_tuple = tuple(sorted(g_clean.items()))
            if g_tuple not in seen:
                seen.add(g_tuple)
                poly_g = sum(c * x**e for e, c in g_clean.items())
                successors.append(Poly(poly_g, x))
                succs_list.append(g_clean)

    return successors, succs_list


def cost(k):
    return math.floor(k + (3 + math.sqrt(8*k - 7)) / 2)


def edge_cost(homotopy_polynomial):
    return sum(v*cost(k) for k, v in homotopy_polynomial.items() if k != 0)


def isInterleaved(chain1, chain2, N):
    """
    Inputs:
        chain1 : list of dict
        chain2 : list of dict
        N : int or float (Number of nodes)

    Outputs:
        bool
            Returns True if the chains are interleaved
            otherwise returns False.
    """

    for i in range(len(chain1) - 1):

        poly_chain1 = Poly(sum(c * x**e for e, c in chain1[i].items()), x)
        poly_chain2 = Poly(sum(c * x**e for e, c in chain2[i].items()), x)


        succs_list1 = succs(chain1[i])[1]
        succs_list2 = succs(chain2[i])[1]


        deriv1 = poly_chain1.diff(x)
        deriv2 = poly_chain2.diff(x)


        m1 = deriv1.eval(1) - poly_chain1.eval(1) + N
        m2 = deriv2.eval(1) - poly_chain2.eval(1) + N


        approved_cost1 = [s for s in succs_list1 if edge_cost(s) <= m1 + 1]
        approved_cost2 = [s for s in succs_list2 if edge_cost(s) <= m2 + 1]


        approved_cost1_tuples = [tuple(sorted(s.items())) for s in approved_cost1]
        approved_cost2_tuples = [tuple(sorted(s.items())) for s in approved_cost2]


        if (tuple(sorted(chain2[i + 1].items())) not in approved_cost1_tuples) and (tuple(sorted(chain1[i + 1].items())) not in approved_cost2_tuples):
            return False

    return True


def homotopy_polynomial(G):
    """
    G: a networkx graph

    Returns:
         (poly, dictinry)
         poly: a sympy.Poly object representation of the homotopy polynomial of G
         dictinry: a dictionary {exponent: coefficient} representation of the homotopy polynomial of G
    """
    cycle_counts = {}
    for idx, component in enumerate(nx.connected_components(G), start=1):
        subG = G.subgraph(component)

        cycles = list(nx.cycle_basis(subG))

        cycles_no_self = [c for c in cycles if len(c) > 1]

        cycle_counts[idx] = len(cycles_no_self)

    poly = Poly(sum(x**v for _, v in cycle_counts.items()), x)

    monoms = poly.monoms()
    coeffs = poly.coeffs()
    dictnry = alpha = dict(zip([m[0] for m in monoms], coeffs))

    return poly, dictnry

def create_chain(list_of_graphs, dictionary=True):
    """
    Inputs:
        list_of_graphs : list of networkx graphs
        dictionary : bool, default True
            - If True, returns dictionary representation of the homotopy polynomial.
            - If False, returns the polynomial object itself.

    Outputs:
        chain : list of dictionaries or sympy.Poly objects
    """
    chain = []
    for i in range(len(list_of_graphs)):
        if dictionary:
            chain.append(homotopy_polynomial(list_of_graphs[i])[1])
        else:
            chain.append(homotopy_polynomial(list_of_graphs[i])[0])

    return chain


def derivative_at_one(poly):
    """Compute derivative of polynomial at x=1."""
    return sum(exp * coeff for exp, coeff in poly.items() if exp != 0)

def evaluate_at(poly, x):
    """Evaluate polynomial at x."""
    return sum(coeff * (x ** exp) for exp, coeff in poly.items())


def homotopy_polynomial_poset(N, succs_func):
    """
    N: initial homotopy polynomial (dict of {exponent: coefficient})
    succs_func: function that takes a polynomial (dict) and returns a list of successors (list of dicts)

    Returns:
        V: list of vertices (homotopy polynomials as dicts)
        E: list of edges (pairs of dicts)
    """
    A = deque([{0: N}])  # Queue of polynomials to process
    V = []  # Vertices
    E = []  # Edges

    while A:
        h = A.pop()  # Select polynomial to process
        if h not in V:  # Not yet processed
            V.append(h)
            # m = dh/dt(1) - h(1) + N
            m = sum(exp * coeff for exp, coeff in h.items() if exp != 0) - sum(coeff for exp, coeff in h.items()) + N #derivative_at_one(h) - evaluate_at(h, 1) + N
            B = succs_func(h)[1]  # Generate successors
            for g in B:
                if edge_cost(g) <= m + 1:  # Edge allowance check
                    A.append(g)
                    E.append((h, g))

    return V, E


def random_edge_chain(N, seed=None):
    """
    Generate a list of graphs starting from N isolated nodes and ending with a complete graph,
    adding one random edge at each step.

    Inputs:
        N : int
            Number of nodes in the graph.
        seed : int or None, optional (default=None)
            Random seed for reproducibility.

    Outputs:
        graph_list : list of networkx.Graph
            A list of graphs where:
            - The first graph has N nodes and no edges.
            - Each subsequent graph has one additional random edge.
            - The last graph is the complete graph on N nodes.
    """
    if seed is not None:
        random.seed(seed)

    # Initialize graph with N nodes and no edges
    G = nx.Graph()
    G.add_nodes_from(range(N))
    graph_list = [G.copy()]

    # All possible edges
    all_edges = [(i, j) for i in range(N) for j in range(i+1, N)]
    random.shuffle(all_edges)

    # Add edges one by one
    for edge in all_edges:
        G.add_edge(*edge)
        graph_list.append(G.copy())

    return graph_list


def label(poly):
    return "{" + ", ".join(f"{k}:{v}" for k,v in sorted(poly.items())) + "}"



def offset_arrow(x0, y0, x1, y1, r_start=0.05, r_end=0.05):
    dx = x1 - x0
    dy = y1 - y0
    dist = math.sqrt(dx ** 2 + dy ** 2)
    if dist == 0:
        return x0, y0, x1, y1  # avoid division by zero

    # unit vector
    ux = dx / dist
    uy = dy / dist

    # offset start and end
    x0_new = x0 + ux * r_start
    y0_new = y0 + uy * r_start
    x1_new = x1 - ux * r_end
    y1_new = y1 - uy * r_end

    return x0_new, y0_new, x1_new, y1_new



def poset_visualization(V, E, r_start=0.42, r_end=0.42, chain1=None, chain2=None):

    G = nx.DiGraph()
    for v in V:
        G.add_node(label(v))
    for u, v in E:
        G.add_edge(label(u), label(v))

    # Convert chain nodes to label sets
    chain1 = set(label(x) for x in chain1) if chain1 else set()
    chain2 = set(label(x) for x in chain2) if chain2 else set()


    sources = [n for n in G.nodes() if G.in_degree(n) == 0]
    source = sources[0] if sources else list(G.nodes())[0]

    levels = {}
    for node in nx.topological_sort(G):
        if node in sources:
            levels[node] = 0
        else:
            levels[node] = 1 + max(levels[p] for p in G.predecessors(node))

    level_groups = defaultdict(list)
    for node, lvl in levels.items():
        level_groups[lvl].append(node)

    pos = {}
    x_spacing = 1.5
    y_spacing = 1.0

    for lvl, nodes_at_level in level_groups.items():
        n_nodes = len(nodes_at_level)
        y_positions = [i * y_spacing for i in range(n_nodes)]
        y_centered = [y - (sum(y_positions) / n_nodes) for y in y_positions]

        for i, node in enumerate(nodes_at_level):
            pos[node] = (lvl * x_spacing, y_centered[i])


    # Normal nodes
    normal_node_trace = go.Scatter(
        x=[pos[k][0] for k in G.nodes()],
        y=[pos[k][1] for k in G.nodes()],
        mode="markers",
        hovertext=list(G.nodes()),
        marker=dict(
            symbol="circle",
            size=20,
            color="red",
            opacity=1.0,
            line=dict(color="#1F1F1F", width=1.2)
        ),
        hoverinfo="text",
        visible=True
    )

    # Full highlighted nodes (both chains)
    full_colors = []
    full_opacities = []

    for n in G.nodes():
        in1 = n in chain1
        in2 = n in chain2
        if in1 and in2:
            color = "#f542a4"     # overlap
        elif in1:
            color = "#0074D9"     # chain1
        elif in2:
            color = "#2ECC40"     # chain2
        else:
            color = "gainsboro"   # not in chains
        full_colors.append(color)
        full_opacities.append(1.0 if (in1 or in2) else 0.75)

    highlight_node_trace = go.Scatter(
        x=[pos[k][0] for k in G.nodes()],
        y=[pos[k][1] for k in G.nodes()],
        mode="markers",
        hovertext=list(G.nodes()),
        marker=dict(
            symbol="circle",
            size=20,
            color=full_colors,
            opacity=full_opacities,
            line=dict(color="#1F1F1F", width=1.2)
        ),
        hoverinfo="text",
        visible=False
    )

    # Only Chain1 nodes
    chain1_colors = ["#0074D9" if n in chain1 else "gainsboro" for n in G.nodes()]
    chain1_opacities = [1.0 if n in chain1 else 0.75 for n in G.nodes()]

    chain1_node_trace = go.Scatter(
        x=[pos[k][0] for k in G.nodes()],
        y=[pos[k][1] for k in G.nodes()],
        mode="markers",
        hovertext=list(G.nodes()),
        marker=dict(
            size=20,
            color=chain1_colors,
            opacity=chain1_opacities,
            line=dict(color="#1F1F1F", width=1.2)
        ),
        hoverinfo="text",
        visible=False
    )

    # Only Chain2 nodes
    chain2_colors = ["#2ECC40" if n in chain2 else "gainsboro" for n in G.nodes()]
    chain2_opacities = [1.0 if n in chain2 else 0.75 for n in G.nodes()]

    chain2_node_trace = go.Scatter(
        x=[pos[k][0] for k in G.nodes()],
        y=[pos[k][1] for k in G.nodes()],
        mode="markers",
        hovertext=list(G.nodes()),
        marker=dict(
            size=20,
            color=chain2_colors,
            opacity=chain2_opacities,
            line=dict(color="#1F1F1F", width=1.2)
        ),
        hoverinfo="text",
        visible=False
    )


    def edge_color(u, v):
        in1 = (u in chain1 and v in chain1)
        in2 = (u in chain2 and v in chain2)
        if in1 and in2:
            return "#f542a4", 1.0
        elif in1:
            return "#0074D9", 1.0
        elif in2:
            return "#2ECC40", 1.0
        return "gainsboro", 0.75

    normal_edges = []
    chain1_edges = []
    chain2_edges = []
    full_edges = []

    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]

        # normal black edges
        normal_edges.append(go.Scatter(
            x=[x0, x1], y=[y0, y1],
            mode="lines",
            line=dict(color="black", width=3),
            hoverinfo="none",
            visible=True
        ))

        # full highlight
        color, op = edge_color(u, v)
        full_edges.append(go.Scatter(
            x=[x0, x1], y=[y0, y1],
            mode="lines",
            line=dict(color=color, width=3),
            opacity=op,
            hoverinfo="none",
            visible=False
        ))

        # only chain1 edges
        if u in chain1 and v in chain1:
            c = "#0074D9"
            o = 1.0
        else:
            c = "gainsboro"
            o = 0.75
        chain1_edges.append(go.Scatter(
            x=[x0, x1], y=[y0, y1],
            mode="lines",
            line=dict(color=c, width=3),
            opacity=o,
            hoverinfo="none",
            visible=False
        ))

        # only chain2 edges
        if u in chain2 and v in chain2:
            c = "#2ECC40"
            o = 1.0
        else:
            c = "gainsboro"
            o = 0.75
        chain2_edges.append(go.Scatter(
            x=[x0, x1], y=[y0, y1],
            mode="lines",
            line=dict(color=c, width=3),
            opacity=o,
            hoverinfo="none",
            visible=False
        ))


    traces = (
        normal_edges + [normal_node_trace] +
        chain1_edges + [chain1_node_trace] +
        chain2_edges + [chain2_node_trace] +
        full_edges + [highlight_node_trace]
    )

    fig = go.Figure(data=traces)

    # trace blocks
    n_normal = len(normal_edges) + 1
    n_c1 = len(chain1_edges) + 1
    n_c2 = len(chain2_edges) + 1
    n_full = len(full_edges) + 1

    # visibility masks
    mask_normal = [True]*n_normal + [False]*(n_c1 + n_c2 + n_full)
    mask_c1     = [False]*n_normal + [True]*n_c1 + [False]*(n_c2 + n_full)
    mask_c2     = [False]*(n_normal + n_c1) + [True]*n_c2 + [False]*n_full
    mask_full   = [False]*(n_normal + n_c1 + n_c2) + [True]*n_full


    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                x=0.01, y=0.99,
                buttons=[
                    dict(label="Normal",         method="update", args=[{"visible": mask_normal}]),
                    dict(label="Only Chain 1",   method="update", args=[{"visible": mask_c1}]),
                    dict(label="Only Chain 2",   method="update", args=[{"visible": mask_c2}]),
                    dict(label="Both Chains",    method="update", args=[{"visible": mask_full}]),
                ]
            )
        ],
        title="Homotopy Polynomial Poset",
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    fig.show()


V, E = homotopy_polynomial_poset(5, succs)

#print(E)
print('Homotopy poset construction...Done')

graphs1 = random_edge_chain(5)
graphs2 = random_edge_chain(5)

#chain1 = [{0: 6}, {0: 5}, {0: 4}, {0: 3}, {0: 2}, {0: 1}, {1: 1}, {2: 1}, {3: 1}, {4: 1}, {5: 1}, {6: 1}, {7: 1}, {8: 1}, {9: 1}, {10: 1}]#create_chain(graphs1, dictionary=True)
#chain2 = [{0: 6}, {0: 5}, {0: 4}, {0: 3}, {1: 1, 0: 2}, {2: 1, 0: 2}, {2: 1, 0: 1}, {2: 1}, {3: 1}, {4: 1}, {5: 1}, {6: 1}, {7: 1}, {8: 1}, {9: 1}, {10: 1}]#create_chain(graphs2, dictionary=True)

chain1 = create_chain(graphs1, dictionary=True)
chain2 = create_chain(graphs2, dictionary=True)
#print(chain1==chain2)

poset_visualization(V, E, r_start=0.42, r_end=0.42, chain1=chain1, chain2=chain2)

print(isInterleaved(chain1, chain2, 5))
