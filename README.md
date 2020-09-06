
# Supergraphs

Implementation of the method described in *Local Patterns and Supergraph for Chemical Graph Classification with Convolutional Networks*
(S+SSPR2018, [pdf here](https://bougleux.users.greyc.fr/articles/sspr18graphconvnet.pdf)).

## List of executables

* bin/supergraph : compute a supergraph of a graph database as well as the projections of input data onto the supergraph
* bin/supergraph-cv : compute a supergraph for each cross-validation subset
* bin/supergraph-cv-omp : the same, multi-threaded with OpenMP
* bin/project : compute the projections of the database's graphs onto a supergraph
* bin/stars : compute stars present in a graph database encoded as a vector
* bin/ds2json : convert a `.ds` dataset into json format
* bin/ds_stats : compute stats of a dataset
* bin/cv_indices : create files of indices for each cross-validation subset


It is possible to create the program `bin/supergraph-cv-omp` using a bipartite approximation of
graph edit distance at each step of the algorithm by using the rule `supergraph-cv-bipartite` of
the Makefile. This accelerates the computation but the resulting supergraph will be more dense.
For more information about this, please consult the paper.

## Requirements

# *lsape*, *graph_lib* and *edmonds*

_lsape is_ a toolbox to address the linear sum assignment problem with edition (Greyc lab).\
_graph_lib_ is toolbox of graph algorithms, especially adressing the problem of Graph Edit Distance estimation (Greyc lab).\
_edmonds_ is an implementation of the Edmonds-Blossom algorithm to find minimum weight perfect matchings on graphs.

```
* lsape                 https://bougleux.users.greyc.fr/lsape
* graph_lib             https://github.com/bgauzere/graph-lib
* edmonds               https://github.com/lacop/edmonds
```

It is then necessary to update the global variables in the Makefile :


* `LSAPE_DIR`           Path to the repertory `include` of *lsape*
* `GRAPH_DIR`           Path to the repertory `include` of *graph_lib*
* `GRAPH_OBJ`           Path to the file `graphlib.a`
* `EDMONDS_DIR`         Path to the root repertory of *edmonds*

for *graph_lib* dependences are :

* LSAPE
* TinyXML
* Eigen3
* OpenMP (for the multi-thread version)

