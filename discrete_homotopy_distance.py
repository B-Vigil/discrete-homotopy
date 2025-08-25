import copy
from sympy import Poly, symbols
from sympy.abc import x
from copy import deepcopy
import math
import networkx as nx

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