from pygsp import graphs
import networkx as nx
import json



def nx_graph( filename ):
    return nx.json_graph.node_link_graph(json.loads(open(filename).read()), directed=False, multigraph=False)

def pg_graph( filename ):
    return graphs.Graph(nx.convert_matrix.to_numpy_matrix( nx_graph(filename) ))

def proj( filename, N ):
    j_proj = json.loads(open( filename ).read())

    p = np.zeros((N, len(j_proj)), dtype=np.int32)
    y = np.array(len(j_proj))

    for i, proj in enumerate(j_proj):
        p[np.array(proj['nodes']), i] = 1
        y[i] = proj['value']

    return p, y


