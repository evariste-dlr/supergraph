import json
import sys, os
sys.path.insert(0, '..')

from lib.dcnn.python import data

A,X,y = data.parse_graph_data()

j_A = []

i=0
for g,n,_y in zip(A,X,y):
    graph = dict()
    graph['graph'] = i
    graph['adj'] = g.astype(int).tolist()
    graph['nodes'] = [x.index(1)+1 for x in n.tolist()]
    graph['value'] = _y.tolist()[0].index(1)
    j_A.append(graph)
    i+=1


with open('nci1.json', 'w') as f:
    json.dump(j_A, f, sort_keys=True)
