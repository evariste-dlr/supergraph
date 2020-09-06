/**
 * @file stars.h
 * @author Évariste Daller <<evariste.daller@unicaen.fr>>
 * @date 2017/12/13
 * @version 0
 * @requirements graph-lib
 *
 * Methods to count 1-hop treelets frequencies over a graph
 * and count the distribution of treelets around each node
 *
 */


#ifndef _TREELETS_H_
#define _TREELETS_H_

#define STAR_CODE_BASE 256  /**< When using a long to encode a star code, base for a node */

#include <vector>
#include <list>
#include <unordered_map>

#include <algorithm>
#include <cstdlib>


#include "graph.h"
#include "Dataset.h"



/*****************************************
 *            Data structures            *
 *****************************************/

#define MAX_NODE_CODE_CHEMO 120     /**< Maximum code number for an atom (119 corresponds to Deuterium) */
#define STAR_CODE_SIZE      100     /**< Size of a star code, ie. length of the string. Related to the maximum degree */
typedef char* star_t;



/*****************************************
 *       Algorithms  and utils           *
 *****************************************/

/**
 * @brief  Recursive procedure to compute the combinations of elements of an array
 * @see combinations(int*, int, std::list<std::vector<int>*>&)
 */
void combinations(int* array, int n, int* C,  int m,  std::list<std::vector<int>*> & L)
{
  for (int i=0; i<n; i++){
    std::vector<int> * newC = new std::vector<int>(m+1);
    for (int j=0; j<m; j++)
      (*newC)[j] = C[j];

    (*newC)[m] = array[i];
    L.push_back(newC);
    combinations(array+i+1, n-i-1, newC->data(), m+1, L);
  }
}


/**
 * @brief  Compute all the possible combinations of the elements in an array
 *
 *  The result is returned in the list L.
 *  No additional memory than the list is used.
 *
 * @param  array        The elements to combine
 * @param  n            Size of the array
 * @param  L            Resulting list of combinations
 */
void combinations(int* array, int n, std::list<std::vector<int>*> & L)
{
  combinations(array, n, NULL, 0, L);
}


void delete_combinations(std::list<std::vector<int>*> & L)
{
  for (std::list<std::vector<int>*>::iterator it=L.begin(); it!=L.end(); it++)
    delete *it;
}


/**
 * @brief A code corresponding to a chemical graph node (an atom)
 */
uint code_chemo( const int& node_value ){
  return node_value;
}


/**
 * @brief A code corresponding to an atom and its incident edge
 */
uint code_chemo( const int& node_value, const int& edge_value ){
  return (edge_value-1) * MAX_NODE_CODE_CHEMO + code_chemo(node_value);
}




/*************************  *********  ****************************/

/**
 * @brief  Unique code of a 1-hop treelet, invariant to neighbors permutations
 *
 *  The returned code is encoded as follow :
 *  * The couple of labels of each neighbor and its incident edge is encoded into a number
 *  * The obtained set is sorted
 *  * The returned code is a number in base <code>base</code> where each digit corresponds to a neighbor and its edge
 *
 * @param  g            The graph to be considered
 * @param  nodes        Nodes indices in the motif around the center. All nodes need to be neighbors to the center
 * @param  center       Index of the center node
 * @param  f_code_label Fonction that returns a different number for each labels of couple (node,edge)
 * @return  A code corresponding to the motif (heap-allocated)
 */
template <class N, class E>
star_t code_star( Graph<N,E> & g,
                  const std::vector<int>& nodes,
                  int center,
                  uint (&f_code_edge)(const N&, const E&),
                  uint (&f_code_node)(const N&) )
{
  int n = nodes.size();
  uint codes[n+1];

  codes[0] = f_code_node(g[center]->attr);
  for (int i=1; i<n+1; i++){
    int j = nodes[i-1];
    codes[i] = f_code_edge(g[j]->attr, g.getEdge(center,j)->attr);
  }

  std::sort(codes+1,  codes+n);

  // Convert code in a C string
  char* code = new char[STAR_CODE_SIZE];
  char* ptr = code;
  for (int i=0; i<n+1; i++){
    int offset = 0;
    if (codes[i] < 10){
      ptr[offset] = '0';
      offset ++;
    }
    if (codes[i] < 100){
      ptr[offset] = '0';
      offset ++;
    }
    sprintf(ptr+offset, "%d", codes[i]);
    ptr+=3;
  }

  return code;
}



/**
 * @brief  Add all the stars centered on the given node to the dictionnary
 *
 *  This method updates the sparse representation of the node <code>center</code>
 *  in the motif space <code>distribution</code>. The keys of this map refer to
 *  values of <code>dict</code>.
 *
 * @param  g            The graph
 * @param  center       The center index of the stars in the graph
 * @param  dict         The dictionnary to update (add the unkown stars and update frequency)
 * @param  distribution A sparse vector representing the distribution of motifs around <code>center</code>
 * @see  code_star
 */
template <class N, class E>
void add_stars( Graph<N,E> & g,
                int center,
                std::unordered_map<std::string, uint> & dict,
                std::unordered_map<uint, uint> & distribution,
                uint (&f_code_edge)(const N&, const E&),
                uint (&f_code_node)(const N&) )
{
  int n = g[center]->Degree();
  int * neighbors = new int[n];

  GEdge<E> *p = g[center]->getIncidentEdges();
  for (int i=0; i<n && p; i++, p=p->Next())
    neighbors[i] = p->IncidentNode();

  std::list<std::vector<int>*> stars;
  combinations( neighbors, n, stars );  // List in stars all the neighbors combinations

  for (std::list<std::vector<int>*>::iterator it=stars.begin(); it != stars.end(); it++){
    star_t rcode = code_star( g, **it, center, f_code_edge, f_code_node );
    std::string code = rcode;

    if (dict.count(code) == 0)
      dict.insert({code, dict.size()});

    if (distribution.count(dict[code]) == 0)
      distribution.insert({dict[code], 1});
    else
      distribution[dict[code]]++;

    delete[] rcode;
  }

  delete_combinations(stars);
  delete[] neighbors;
}



/**
 * @brief  Compute the dictionnary of the present stars in the given graph
 *
 * @param  g            The graph
 * @param  dict         Dictionnary populated by this function
 * @param  distribution The distribution of motifs over the nodes, in sparse representation
 * @param  f_code_label Function that returns a unique code for a given configuration of labels
 * @see add_stars
 * @see code_star
 */
template <class N, class E>
std::unordered_map<std::string, uint> &
star_dictionnary( Graph<N,E> & g,
                  std::unordered_map<std::string, uint> & dict,
                  std::vector<std::unordered_map<uint,uint> > & distribution,
                  uint (&f_code_edge)(const N&, const E&),
                  uint (&f_code_node)(const N&) )
{
  //if (double(sizeof(long)) / (std::log(base)/std::log(2)) > double(g.Degree()))
  //  throw std::exception("The code of a graph is bigger than a long");
  //TODO omp parallel_for
  for (int i=0; i<g.Size(); i++)
    add_stars(g, i, dict, distribution[i], f_code_edge, f_code_node);
  return dict;
}



/**
 * @brief  Compute and return the dictionnary of the present stars in the given dataset
 *
 * @param   dataset             The database
 * @param   distribution        An array of distribution over the nodes of each graph. Must be allocated before
 * @param   f_code_edge         A reference to a function that computes a code for an edge and the incident node
 * @param   f_code_node         A reference to a function that computes a code for a node
 * @return  The dictionnary of motifs, as a map allocated on the heap
 */
template <class N, class E, typename D>
std::unordered_map<std::string, uint> *
star_dictionnary( Dataset<N,E,D> & dataset,
                  std::vector<std::unordered_map<uint,uint> > distribution[],
                  uint (&f_code_edge)(const N&, const E&),
                  uint (&f_code_node)(const N&) )
{
  std::unordered_map<std::string, uint> * dict = new std::unordered_map<std::string, uint>();
  for (int g=0; g<dataset.size(); g++)
    star_dictionnary(*(dataset[g]), *dict, distribution[g], f_code_edge, f_code_node);

  return dict;
}



/**
 * @brief  Compute and return the dictionnary of stars and their distribution in a Chemical dataset
 */
std::unordered_map<std::string, uint> *
star_dictionnary( Dataset<int,int,double> & dataset,
                  std::vector<std::unordered_map<uint,uint> > distribution[] )
{
  return star_dictionnary( dataset, distribution, code_chemo, code_chemo );
}



#endif
