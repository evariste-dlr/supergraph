/**
 * @file supergraph.h
 * @author Évariste Daller <<evariste.daller@unicaen.fr>>
 * @date 2017/09/08
 * @version 0
 *
 * Set of classes and functions to manage the creation of a supergraph
 * of several graphs. The main purpose is to generate the network's
 * input graph
 *
 * All graphs G<N,E> = (Vertex, Edges) are labeled with nodes attributes
 * of type N and Edges attributes of type E
 *
 * If macro NESTED_TYPE is defined, the behaviour of OpenMP parallelism
 * will be controlled as follow :
 * * NESTED_TYPE=0  :  equiv. NESTED_TYPE not defined, multi-threading for GED computation (mIPFP)
 * * NESTED_TYPE=1  :  Nested parallelism : computation of distance matrices in parallel and GED parallel too
 * * NESTED_TYPE=2  :  Distances in parallel only, the mIPFP will be sequential
 *
 */

#ifndef __SUPERGRAPH_H__
#define __SUPERGRAPH_H__

#ifdef VERBOSE
 #include <time.h>  // For printing time and/or progress bars
#endif

#ifndef NESTED_TYPE
 #define NESTED_TYPE 0
#endif

#if NESTED_TYPE != 0
 #include <omp.h>

 #ifndef NUM_THREADS
  #define NUM_THREADS 16
 #endif
#endif

#include <set>
#include <limits>

// From Graph-lib
#include "graph.h"
#include "GraphEditDistance.h"
#include "Dataset.h"

// From Edmond's algorithm directory
#include "EdmondsMatching.hpp"

#include "SupergraphCostFunction.h"
#include "utils.h"
#include "ECMapping.h"


//#undef VERBOSE




/**
 * Construct a supergraph of g1 and g2 based on the edit path
 * returned by ged(g1,g2). The rule to labelize substituted
 * elements i and j between g1 and g2 in this edit path is given by
 * f_compnodes(g1.node(i),g2.node(j)) and f_compedges(g1.edge(i),g2.edge(j))
 *
 * @brief Construct a supergraph of 2 given graphs
 * @return a new graph containing g1 and g2
 */
template <class N, class E>
Graph<N,E>* supergraph(Graph<N,E> & g1,
                       Graph<N,E> & g2,
                       N (&f_compnodes)(const N&, const N&),
                       E (&f_compedges)(const E&, const E&),
                       GraphEditDistance<N,E> & ged
                       )
{
  int n=g1.Size();
  int m=g2.Size();
  int * rho = new int[n+1];
  int * varrho = new int[m+1];

  ged.getOptimalMapping(&g1, &g2, rho, varrho);
  Graph<N,E>* result = supergraph(g1, g2, f_compnodes, f_compedges, rho, varrho);

  delete[] rho;
  delete[] varrho;

  return result;
}




/**
 * @brief Construct a supergraph of g1 and g2 based on the given error-correcting mapping
 *
 *   The supergraph is constructed as follow :
 *   * Construct all elements in g1 that are not deleted in the mapping
 *   * Construct all insterted ellements in g2
 *   The labels of the common elements between g1 and g2 are deduced from
 *   the functions <code>f_compnodes</code> and <code>f_compedges</code>
 *
 *  @param g1  The first graph
 *  @param g2  The second graph
 *  @param f_compnodes  A pointer to a function comparing a node label from g1 to a node label from g2
 *  @param f_compedges  A pointer to a function comparing an edge label from g1 to another from g2
 *  @return A supergraph of g1 and g2
 */
template <class N, class E>
Graph<N,E>* supergraph(Graph<N,E> & g1,
                       Graph<N,E> & g2,
                       N (&f_compnodes)(const N&, const N&),
                       E (&f_compedges)(const E&, const E&),
                       int rho[], int varrho[])
{
  int n=g1.Size();
  int m=g2.Size();

  Graph<N,E>* sg = new Graph<N,E>(g1.isDirected() && g2.isDirected());


  // Insert all ellements of g1 :
  // i) nodes
  for (int i=0; i<n; i++){
    if (rho[i] < m)
      sg->Add(new GNode<N,E>(i, f_compnodes(g1[i]->attr, g2[rho[i]]->attr)));
    else
      sg->Add(new GNode<N,E>(i, g1[i]->attr));
  }

  // ii) edges
  for (int i=0; i<n; i++){
    GEdge<E> *p = g1[i]->getIncidentEdges();
    while(p){
      int start = i;
      int end = p->IncidentNode();
      int f_start = rho[start];
      int f_end = rho[end];

      if(!(sg->isLinked(start, end))){
        if(( f_start < m) && (f_end < m)){
          GEdge<E>  *mappedEdge = g2.getEdge(f_start,f_end);
          if( mappedEdge != NULL){
            // Substitution -> insert edge with comparison between g1 and g2
            sg->Link(start, end, f_compedges(p->attr, mappedEdge->attr));
          }else{
            //Edge does not exist in G2 : Deletion -> insert edge as in g1
            sg->Link(start, end, p->attr);
          }
        }else{
          //start or end has been deleted : Deletion -> insert as in g1
          sg->Link(start, end, p->attr);
        }
      }
      p = p->Next();
    }
  }


  // Insert elements of g2 :
  //  nodes inserted and their incident edges

  // i) nodes
  int nbNodes = n;
  std::map<int,int> g2nodes;
  for (int j=0; j<m; j++){
    if (varrho[j] >= n){
      sg->Add(new GNode<N,E>(nbNodes, g2[j]->attr));
      g2nodes[j] = nbNodes;
      nbNodes++;
    }
  }

  // ii) edges from common nodes
  for (int j=0; j<m; j++){
    if (varrho[j] < n){
      GEdge<E> *p = g2[j]->getIncidentEdges();
      while(p){
        int g2_start = j;
        int g2_end = p->IncidentNode();
        int g1_start = varrho[g2_start];
        int g1_end = varrho[g2_end];
        if (g1_end >= n && !(sg->isLinked(g1_start, g2nodes[g2_end])))
          sg->Link(g1_start, g2nodes[g2_end], p->attr);
        else if(!(g1.isLinked(g1_start, g1_end)) && !(sg->isLinked(g1_start, g1_end)))
          sg->Link(g1_start, g1_end, p->attr);
        p = p->Next();
      }
    }
  }

  // iii) edges from inserted nodes into g2
  for (std::map<int,int>::iterator it=g2nodes.begin(); it!=g2nodes.end(); it++){
    GEdge<E> *p = g2[it->first]->getIncidentEdges();
    while(p){
      int g2_end = p->IncidentNode();
      int g1_start = it->second;
      int g1_end = varrho[g2_end];
      if(g1_end < n && !(sg->isLinked(g1_start,g1_end))){
        sg->Link(it->second, g1_end, p->attr);
      }
      else if (g1_end >= n && !(sg->isLinked(g1_start, g2nodes[g2_end]))){
        sg->Link(it->second, g2nodes[g2_end], p->attr);
      }
      p = p->Next();
    }
  }

  return sg;
}



/**
 * @brief Construct a supergraph of an entire dataset
 *
 *   Description de la méthode
 *
 * @param  dataset  A reference to the database of graphs, il will be destroyed by the procedure and contains only the returned graph
 * @param f_compnodes  A pointer to a function comparing a node label from g1 to a node label from g2
 * @param f_compedges  A pointer to a function comparing an edge label from g1 to another from g2
 * @param ged  The Graph Edit Distance implementation used to compute the mapping
 * @param ksol   Number of nested threads for ged computation. If There is no parallelism for GED computation, this argument is ignored
 * @param use_bipartite_matrix    If true, a bipartite approximation of the GED is used instead of *ged* to compute each distance matrix
 * @return A pointer to the supergraph, that is, the only remaining graph in the dataset.
 * @note As the supergraph is referenced by the dataset, it will be destroyed in the dataset's destructor.
 */
template <class N, class E>
Graph<N,E>* supergraph( Dataset<N,E,double> & dataset,
                        N (&f_compnodes)(const N&, const N&),
                        E (&f_compedges)(const E&, const E&),
                        GraphEditDistance<N,E> & ged,
                        unsigned int ksol = 10,
                        bool use_bipartite_matrix = false )
{
  int n = dataset.size();
  int count = 0;
  Graph<N,E>* gs = NULL; // surviving graph if n is odd

  /* Use bipartite ged matrix ? */
  GraphEditDistance<N,E>* matrixGED = &ged;
  if (use_bipartite_matrix){
    matrixGED = new BipartiteGraphEditDistanceMulti<N,E> (ged.getCostFunction(), ksol);
  }

  while (n > 1 || gs != NULL){

    if (n%2 != 0){
      if (gs == NULL){ // gs is empty => save the last graph in it
        gs = dataset.pop_back();
        n--;
      }
      else { // gs is already a reserved graph => use it
        dataset.add(gs, 0);
        gs = NULL;
        n++;
      }
    }

    count++;
    std::cout << std::endl;
    std::cout << "  Step " << count << "  (" << n << " graphs)" << std::endl;
    std::cout << "  Distances computation..." << std::endl;


    double elapsed;
    struct timeval t0, t1;

    #ifdef VERBOSE
    time_t tv, tvprev;
    time(&tvprev);
    gettimeofday(&t0, NULL);
    #endif

    // Distance matrix : this is also the adjacency matrix of the graph in which we want to find a mapping
    double * D = new double[n*n];

    #if NESTED_TYPE == 1
     omp_set_nested(1);
    #endif
    #if NESTED_TYPE != 0 || defined BIPARTITE_MATRIX
     int num_group_threads = NUM_THREADS/ksol;
     if (use_bipartite_matrix)  num_group_threads = NUM_THREADS;
     //omp_set_num_threads(NUM_THREADS);
     std::vector<GraphEditDistance<N,E> *> local_ged(num_group_threads);
     for (unsigned int _g=0; _g<local_ged.size(); _g++)
       local_ged[_g] = matrixGED->clone();

     unsigned int completeSteps = 0;
     unsigned int NSteps = (n*(n+1))/2;

     #pragma omp parallel num_threads(num_group_threads)
     {
     #pragma omp for schedule(dynamic) nowait
      for(int k=0; k<n*(n+1)/2; k++) {
          int i = k/(n+1), j = k%(n+1);
          if(j>i) i = n - i - 1, j = n - j;

          #pragma omp atomic
            completeSteps++;

          if (i != j){
            GraphEditDistance<N,E> * lged = local_ged[omp_get_thread_num()];
            D[sub2ind(j,i,n)] = D[sub2ind(i,j,n)] = std::min((*lged)(dataset[i], dataset[j]), (*lged)(dataset[j], dataset[i]));

           #ifdef VERBOSE
            if(omp_get_thread_num() == 0)
            {
              time(&tv);
              if (difftime(tv, tvprev) >= 1.0){
                time(&tvprev);
                int progress = 51*(float)(completeSteps)/NSteps;
                std::cout << "\r" << "[" << std::string(progress, '=') << std::string(std::max(0,50-progress), ' ') << "]  (" << 2*progress << "%)";
                std::cout.flush();
              }
            }
           #endif
          }
      }
     }// end omp parallel

     for  (unsigned int _g=0; _g<local_ged.size(); _g++)
       delete local_ged[_g];

    #else
    for (int j=0; j<n; j++){
      for (int i=0; i<j; i++){
        D[sub2ind(j,i,n)] = D[sub2ind(i,j,n)] = std::min((*matrixGED)(dataset[i], dataset[j]), (*matrixGED)(dataset[j], dataset[i]));
        #ifdef VERBOSE
        time(&tv);
        if (difftime(tv, tvprev) >= 1.0){
          time(&tvprev);
          int progress = 51*(float)(j*n+i)/(n*n);
          std::cout << "\r" << "[" << std::string(progress, '=') << std::string(50-progress, ' ') << "]  (" << 2*progress << "%)";
          std::cout.flush();
        }
        #endif
      }
    }
    #endif // NESTED_TYPE != 0
    #ifdef VERBOSE
    gettimeofday(&t1, NULL);
    elapsed = (double)(t1.tv_usec - t0.tv_usec)/1000000 + (double)(t1.tv_sec - t0.tv_sec);
    std::cout << "\r  [" << GREEN << "done (" << elapsed << " s)" << RESET << "]" << std::string(55, ' ') << std::endl;
    #endif

    for (int i=0; i<n; i++)
      D[sub2ind(i,i,n)] = std::numeric_limits<double>::max();

    EdmondsMatching edmond(n);
    edmond.fromMatrix(D, n);
    edmond.findMinimumWeightMatching();
    auto edges = edmond.getMatchedEdges();

    // Generate n/2 supergraphs
    Graph<N,E>* garray[n/2];
    int i=0;
    std::cout << "  Computing supergraphs..." << std::endl;
    gettimeofday(&t0, NULL);
    for (auto e: edges){
      //std::cerr << "  (" << e.first << " ; " << e.second << ")" << std::endl;
      garray[i] = supergraph(*(dataset[e.first]), *(dataset[e.second]),
                             f_compnodes, f_compedges, ged);
      i++;
    }
    gettimeofday(&t1, NULL);
    elapsed = (double)(t1.tv_usec - t0.tv_usec)/1000000 + (double)(t1.tv_sec - t0.tv_sec);
    std::cout << "\r  [" << GREEN << "done (" << elapsed << " s)" << RESET << "]" << std::endl;

    // Update the dataset with the new graphs
    dataset.clear();
    for (i=0; i<n/2; i++){
      dataset.add(garray[i], 0);
    }

    n = n/2;
    delete [] D;

  } // end while

  if (use_bipartite_matrix)
    delete matrixGED;

  return dataset[0];
}



template<class N, class E>
double projectionCost( Graph<N,E>& g,
                       Graph<N,E>& s,
                       int* rho, int* varrho,
                       GraphEditDistance<N,E>& ged,
                       EditDistanceCost<N,E>* proj_cf )
{
  EditDistanceCost<N,E> * pcf;
  EditDistanceCost<N,E> * cf = ged.getCostFunction();
  if (proj_cf == NULL)
    pcf = new ProjectionDistanceCost<N,E>(cf);
  else
    pcf = proj_cf;

  ged.setCostFunction(pcf);
  double dist = ged.GedFromMapping(&g, &s, rho, g.Size(), varrho, s.Size());
  ged.setCostFunction(cf);

  if (proj_cf == NULL) delete pcf;
  return dist;
}



/**
 *  @brief A projection on the supergraph
 *
 *    The projection is realized trough the computation of
 *    a mapping from \ref g to \ref s according to the
 *    given graph edit distance method \ref ged
 *
 * @param g   A graph to be projected
 * @param s   The supergraph
 * @param ged The distance to compute a mapping
 * @param mapping  An array of the size of \reg g contaning the returned projection (where are projected the nodes of g
 * @param proj_cf  The edit cost function to use to compute the cost of the resulting mapping. NULL to use a standard ProjectionDistanceCost with zero insertion and deletion costs
 * @return the corresponding cost (substitution cost)
 * @note \ref mapping must be allocated first
 */
template<class N, class E>
double projection( Graph<N,E> & g,
                   Graph<N,E> & s,
                   GraphEditDistance<N,E> & ged,
                   int * mapping,
                   EditDistanceCost<N,E> * proj_cf = NULL )
{
  int n = g.Size();
  int m = s.Size();
  int* rho = new int[n+1];
  int* varrho = new int[m+1];

  ged.getOptimalMapping(&g, &s, rho, varrho);

  for (int i=0; i<n; i++){
    mapping[i] = rho[i];
  }

  double dist = projectionCost<N,E>(g, s, rho, varrho, ged, proj_cf);
  delete [] rho;
  delete [] varrho;

  return dist;
}


/**
 *  @brief A set of projections based on multiple mappings of the same approximate GED cost
 *
 *    The mappings used as projections are the refined mappings returned by \ref ged.
 *    Identical mappings are then removed avoiding identical projections in the result of
 *    this function.
 *
 * @param g   A graph to be projected
 * @param s   The supergraph
 * @param ged The distance to compute a mapping
 * @param mapping  A list of arrays, each of the size of \reg g contaning the returned projections (where are projected the nodes of g
 * @param proj_cf  The edit cost function to use to compute the cost of the resulting mapping. NULL to use a standard ProjectionDistanceCost with zero insertion and deletion costs
 * @param data_augment  Percentage of data augmentation by considering more sub-optimal mappings
 * @return An array of size mappings.size(), containning the corresponding costs (substitution cost)
 * @note \ref mappings are allocated internally on the heap, and memory management is left to the user
 */
template<class N, class E>
double  projection( Graph<N,E> & g,
                    Graph<N,E> & s,
                    MultistartRefinementGraphEditDistance<N,E> & ged,
                    std::list<ECMapping*> & mappings,
                    EditDistanceCost<N,E> * proj_cf = NULL,
                    double data_augment = .0 )
{
  int n = g.Size();
  int m = s.Size();
  
  std::list<int*> _refinedMappings = ged.getBetterMappings(&g, &s);
  std::list<int*> _reverseMappings = ged.getReverseMappings();

  for ( std::list<int*>::iterator it = _refinedMappings.begin(), it2 = _reverseMappings.begin();
        it != _refinedMappings.end() && it != _reverseMappings.end();
        it++, it2++ ){
    ECMapping* map = new ECMapping(*it, *it2, n, m, ged.GedFromMapping(&g, &s, *it, n, *it2, m));
    mappings.push_back(map);
  }


  // Remove identical mappings
  for( std::list<ECMapping*>::iterator it = mappings.begin();
       it != mappings.end(); it++ ){
    std::list<ECMapping*>::iterator jt = std::next(it);

    while(jt != mappings.end()){
      bool eq = true;

      for (int i=0; eq && i<n; i++)
        if ((*it)->f(i) != (*jt)->f(i))
          eq = false;

      for (int i=0; eq && i<m; i++)
        if ((*it)->r(i) != (*jt)->r(i))
          eq = false;

      if (eq){
        delete (*jt);
        jt = mappings.erase(jt);
      }
      else
        jt++;
    }
  }


  // Compute min cost
  double minCost=std::numeric_limits<double>::max();
  for (std::list<ECMapping*>::const_iterator it = mappings.begin(); it != mappings.end(); it++)
    if (minCost > (*it)->cost()) minCost = (*it)->cost();


  // -> keep only best mappings + data augmentation
  mappings.sort(ECMapping::leq_ptr);
  std::list<ECMapping*>::iterator it = mappings.begin();
  bool augmented=false;
  int countmax = mappings.size();
  int i=0;
  while (it != mappings.end()){
    if ((*it)->cost() > minCost && !augmented){
      augmented = true;
      countmax = std::min((double)(i) + data_augment*i, (double)(countmax));
    }
    if (i >= countmax){
      delete (*it);
      it = mappings.erase(it);
    }
    else it++;
    i++;
  }

  return minCost;
}

/****************************************
               Utilitaires
 ****************************************/


std::set<int>
compNode_union(const std::set<int> & first,  const std::set<int> & sec)
{
  std::set<int> output(first);
  for(std::set<int>::iterator it=sec.begin(); it!=sec.end(); it++){
    output.insert(*it);
  }
  return output;
}


double
compNode_mean(const double& first, const double& sec){
  return (first + sec) / 2;
}


template <class N>
N compNode_first(const N& first, const N& sec){
    return sec;
}

template <class N>
N compNode_carbon(const N& first, const N& sec){
    if (first == sec && first == 6)
      return 0;
    else return 1;
}

template <class E>
E compEdge_first(const E& first, const E& sec){
    return first;
}


#endif //__SUPERGRAPH_H__
