import copy
from sympy import Poly, symbols
from sympy.abc import x
from copy import deepcopy
import math
import networkx as nx

import plotly.graph_objects as go

from collections import deque

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
        chain1 : list of dict
        chain2 : list of dict
        N : int or float (Number of nodes)

    Returns:
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
        list_of_graphs : list of networkx graphs
        dictionary : bool, default True
            - If True, returns dictionary representation of the homotopy polynomial.
            - If False, returns the polynomial object itself.

    Returns:
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
    A = deque([{0: N}])  
    V = []  
    E = []  
    while A:
        h = A.pop()  
        if h not in V:  
            V.append(h)
            # m = dh/dt(1) - h(1) + N
            m = sum(exp * coeff for exp, coeff in h.items() if exp != 0) - sum(coeff for exp, coeff in h.items()) + N 
            B = succs_func(h)[1] 
            for g in B:
                if edge_cost(g) <= m + 1: 
                    A.append(g)
                    E.append((h, g))


    return V, E



def offset_arrow(x0, y0, x1, y1, r_start=0.05, r_end=0.05):
    dx = x1- x0
    dy = y1 - y0

    dist = math.sqrt(dx ** 2 + dy ** 2)
    if dist == 0:
        return x0, y0, x1, y1
    
    ux = dx / dist
    uy = dy / dist

    x0_new = x0 + ux * r_start
    y0_new = y0 + uy * r_start
    x1_new = x1 + ux * r_end
    y1_new = y1 + uy * r_end

    return x0_new, y0_new, x1_new, y1_new


def label(poly):
    return "{" + ", ".join(f"{k}:{v}" for k,v in sorted(poly.items())) + "}"


def poset_visualization(V, E, r_start=0.05, r_end=0.05):
    
    G = nx.DiGraph()
    for v in V:
        G.add_node(label(v))
    for u, v in E:
        G.add_edge(label(u), label(v))

    sources = [n for n in G.nodes() if G.in_degree(n) == 0]
    source = sources[0] if sources else list(G.nodes())[0]

   
    levels = {}
    for node in nx.topological_sort(G):
        if node in sources:
            levels[node] = 0
        else:
            levels[node] = 1 + max(levels[pred] for pred in G.predecessors(node))

    
    from collections import defaultdict
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

    
    x_nodes = [pos[k][0] for k in G.nodes()]
    y_nodes = [pos[k][1] for k in G.nodes()]

    node_trace = go.Scatter(
        x=x_nodes, y=y_nodes,
        mode='markers',
        #mode='markers+text',
        #text=list(G.nodes()),
        #textposition='top center',
        hovertext=list(G.nodes()),
        marker=dict(
            symbol='circle',
            size=20,
            color='#FD3216',
            line=dict(color='#1F1F1F', width=1.2)
        ),
        hoverinfo='text'
    )

    
    edge_traces = []
    annotations = []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_traces.append(go.Scatter(
            x=[x0, x1],
            y=[y0, y1],
            mode='lines',
            line=dict(color='black', width=1.5),
            hoverinfo='none'
        ))
        x0_adj, y0_adj, x1_adj, y1_adj = offset_arrow(x0, y0, x1, y1, r_start, r_end)
        annotations.append(dict(
            ax=x0_adj, ay=y0_adj, x=x1_adj, y=y1_adj,
            xref='x', yref='y', axref='x', ayref='y',
            showarrow=True, arrowhead=3, arrowsize=1.8,
            arrowwidth=1.4, arrowcolor='black', opacity=0.7
        ))

    fig = go.Figure(data=edge_traces + [node_trace])
    fig.update_layout(
        title='Layered (Columnar) Graph Layout',
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=0, r=0, b=0, t=40),
        annotations=annotations
    )

    fig.show()
