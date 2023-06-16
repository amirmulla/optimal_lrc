from sage.all import *
from networkx.algorithms.connectivity import EdgeComponentAuxGraph

def get_error_position_subgraph(G, h_idx=0):
    aux_graph = EdgeComponentAuxGraph.construct(G)
    l_subgraphs = sorted(map(sorted, aux_graph.k_edge_components(k=1)))
    tmp = []
    for n in l_subgraphs[h_idx]:
        if G.nodes[n]['bipartite'] == 0:
            tmp.append(G.nodes[n]['coset'])

    l_evalpts_idx = flatten(tmp)
    l_evalpts_idx = [int(str(x)) for x in l_evalpts_idx]
    l_evalpts_idx.sort()

    return l_evalpts_idx