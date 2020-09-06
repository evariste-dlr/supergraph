
#ifndef __SUPERGRAPHCOSTFUNCTION_H__
#define __SUPERGRAPHCOSTFUNCTION_H__

#include <set>
#include "GraphEditDistance.h"


/**
 * @brief A \f$\Chi^2$ distance between sets of integer
 *
 *  The substitution costs between 2 sets \f$S_1\text{ and } S_2$ are computed as :
 *  \f[
 *    d^2(S_1, S_2)  = \sum_{e\in S_1 \cup S_2} \text{freq}(e) \Big(\frac{S_1^e}{|S_1|} - \frac{S_2^e}{|S_2|}\Big)^2
 *  \f]
 *  where freq\f$(e)$ is the proportion of the number of occurences of \f$e$ in all the sets in the two graphs, 
 *  and \f$S^e = 1$ if and only if \f$e\in S$, else it is 0.
 */
class ChiSquareDistanceCost : public EditDistanceCost<std::set<int>, int>
{

 protected:
  double cni;
  double cnd;
  double cei;
  double ced;

  // We store pointers to the graphs on which we work, to know if the histogram is up to date
  Graph<std::set<int>,int> * g1; 
  Graph<std::set<int>,int> * g2;

  std::map<int, int> histogram; //!< Store the number of each label in the 2 graphs
  int Nlabels; //!< Total number of labels in all the sets (in the 2 graphs)


 protected:

  /**
   * @brief Check if the graphs we are working on have changed, if so, re-compute the histogram
   */
  void checkGraphChanged(Graph<std::set<int>,int> * g1, Graph<std::set<int>,int> * g2){
    if (g1 != this->g1 || g2 != this->g2){
      histogram.clear();
      Nlabels = 0;
      computeHistogram(g1, g2);
    }
  }


  /**
   * @brief Compute a the label distribution in the graph(s) as a histogram
   * @see histogram
   */
  void computeHistogram(Graph<std::set<int>,int> * g1, Graph<std::set<int>,int> * g2 = NULL){
    int n = g1->Size();
    int m = (g2 == NULL) ? 0 : g2->Size();
    Nlabels = 0;

    for (int i=0; i<n; i++){
      std::set<int> * set = &((*g1)[i]->attr);
      for (std::set<int>::iterator it=set->begin(); it!=set->end(); it++){
	try{ histogram.at(*it) += 1; }
	catch(std::exception& e){
	  histogram[*it] = 1;
	}
	Nlabels ++;
      }
    }

    for (int i=0; i<m; i++){
      std::set<int> * set = &((*g2)[i]->attr);
      for (std::set<int>::iterator it=set->begin(); it!=set->end(); it++){
	try{ histogram.at(*it) += 1; }
	catch(std::exception& e){
	  histogram[*it] = 1;
	}
	Nlabels ++;
      }
    }
  }

  

 public:

  /**
   * @brief Returns the Chi square distance between n1 and n2
   */
  virtual double NodeSubstitutionCost( GNode<std::set<int>,int> * n1,
				       GNode<std::set<int>,int> * n2,
				       Graph<std::set<int>,int> * g1,
				       Graph<std::set<int>,int> * g2 )
  {
    if (n1->attr == n2->attr)
      return 0;
    else {
      checkGraphChanged(g1, g2);
      
      double dist = 0;
      std::set<int> * set1 = &(n1->attr);
      std::set<int> * set2 = &(n2->attr);
      std::set<int> setUnion(*set1);
            
      for (std::set<int>::iterator it=set2->begin(); it!=set2->end(); it++)
	setUnion.insert(*it);
      

      for (std::set<int>::iterator it=setUnion.begin(); it!=setUnion.end(); it++){
	try{
	  double freq = (double) (histogram.at(*it)) / Nlabels;
	  double term = .0;
	  if (set1->count(*it) != 0) term =  1/set1->size();
	  if (set2->count(*it) != 0) term -= 1/set2->size();
	  term *= term;
	  dist += term * freq;
	}
	catch(std::exception& e){
	  std::cerr << "[E] in NodeSubstitutionCost : histogram does not contain " << *it << ", wtf ?" << std::endl;
	}
      }
      return dist;
    }
  }
  
  virtual double NodeDeletionCost(GNode<std::set<int>,int> * n1,Graph<std::set<int>,int> * g1){
    return cnd;
  }
  
  virtual double NodeInsertionCost(GNode<std::set<int>,int> * n2,Graph<std::set<int>,int> * g2){
    return cni;
  }
  
  
  
  virtual double EdgeSubstitutionCost(GEdge<int> * e1,GEdge<int> * e2,Graph<std::set<int>,int> * g1,Graph<std::set<int>,int> * g2){
    return 1;
  }
  
  virtual double EdgeDeletionCost(GEdge<int> * e1,Graph<std::set<int>,int> * g1){
    return ced;
  }
  
  virtual double EdgeInsertionCost(GEdge<int> * e2,Graph<std::set<int>,int> * g2){
    return cei;
  }

 public:

  /**
   * @brief Construct a cost function where insertion and deletion costs are 1
   */
  ChiSquareDistanceCost () :
    cni(1.0), cnd(1.0),
    cei(1.0), ced(1.0),
    g1(NULL), g2(NULL), Nlabels(0)
  {}

  ChiSquareDistanceCost (double ni,  double nd,  double ei,  double ed) :
    cni(ni), cnd(nd),
    cei(ei), ced(ed),
    g1(NULL), g2(NULL), Nlabels(0)
  {}


 public:

    virtual ChiSquareDistanceCost * clone() const { return new ChiSquareDistanceCost(*this); }

};





/**
 * @brief Defines a projection cost of a graph on another
 */
template<class N, class E>
  class ProjectionDistanceCost : public EditDistanceCost<N,E>
{
 protected:
  EditDistanceCost<N,E> * cf_ref;

  double cni;
  double cnd;
  double cei;
  double ced;

 public:

  virtual double NodeSubstitutionCost(GNode<N,E> * n1, GNode<N,E> * n2, Graph<N,E> * g1, Graph<N,E> * g2){
    return cf_ref->NodeSubstitutionCost(n1, n2, g1, g2);
  }
  
  virtual double NodeDeletionCost(GNode<N,E> * n1, Graph<N,E> * g1){ return cnd; }
  
  virtual double NodeInsertionCost(GNode<N,E> * n2,Graph<N,E> * g2){ return cni; }
  
  
  
  virtual double EdgeSubstitutionCost(GEdge<int> * e1,GEdge<int> * e2,Graph<N,E> * g1,Graph<N,E> * g2){
    return cf_ref->EdgeSubstitutionCost(e1, e2, g1, g2);
  }
  
  virtual double EdgeDeletionCost(GEdge<int> * e1,Graph<N,E> * g1){ return ced; }
  
  virtual double EdgeInsertionCost(GEdge<int> * e2,Graph<N,E> * g2){ return cei; }

  
  
 public:

  /**
   * @brief A projection edit distance cost with zero-valued insertion and deletion costs
   * 
   *   The projection cost of an edit distance mapping is the sum
   *   of the costs of substituted elements. Thus, the substitution cost
   *   of the given existing cost function will be used.
   *   The deletion and insertion costs will be 0.
   */
  ProjectionDistanceCost(EditDistanceCost<N,E> * cf):
    cf_ref(cf),
    cni(.0), cnd(.0),
    cei(.0), ced(.0)
  {}

  /**
   * @brief Initialize a projection edit distance cost with non-zero insertion and deletion costs
   *
   * @param cf  A pointer to the cost function to use for substitution costs
   * @param ni  Cost of node insertions
   * @param nd  Cost of node deletions
   * @param ei  Cost of edge insertions
   * @param ed  Cost of edge deletions
   */
  ProjectionDistanceCost(EditDistanceCost<N,E> * cf, double ni, double nd, double ei, double ed):
    cf_ref(cf),
    cni(ni), cnd(nd),
    cei(ei), ced(ed)
  {}

  virtual ProjectionDistanceCost<N,E> * clone() const { return new ProjectionDistanceCost(*this); }

};




/**
 * @brief A Cost function that allows to construct a Common Supergraph of 2 graphs
 *        without any substitutions, ie. with no distortion in terms of labels
 *
 *   
 */
template <class N, class E>
class CSCostFunction : public EditDistanceCost<N,E>{

protected:
  Graph<N,E> * g1;
  Graph<N,E> * g2;
  int max_degree;

  
protected:
  int maxDegree(Graph<N,E> * g){
    int md=0;
    int d;
    for (int i = 0; i<g->Size(); i++){
      d = (*g)[i]->Degree();
      if (max_degree < d) max_degree = d;
    }
    return md;
  }

  bool checkG1(Graph<N,E> * g){
    if (g1 != g){
      g1 = g;
      max_degree = maxDegree(g1);
      return true;
    }
    return false;
  }

  bool checkG2(Graph<N,E> * g){
    if (g2 != g){
      g2 = g;
      max_degree = maxDegree(g2);
      return true;
    }
    return false;
  }

  bool checkGraphs(Graph<N,E> * g1, Graph<N,E> * g2){
    return checkG1(g1) | checkG2(g2);
  }

  
public:

  CSCostFunction():
    g1(NULL),
    g2(NULL),
    max_degree(0)
  {}

  
  CSCostFunction( const CSCostFunction& other ):
    g1(other.g1),
    g2(other.g2),
    max_degree(other.max_degree)
  {}
  

  virtual EditDistanceCost<N,E> * clone() const {
    return new CSCostFunction(*this);
  }
  
public:

  /**
   * Cost of non-equal node substitution :
   * \f$ 3 + 2\text{max\_degree}\f$ <em>ie.</em> Cost of deleting the node, re-insert it, delete all
   * incident egdes in the worst case and re-insert all of it + 1
   */
  virtual double NodeSubstitutionCost( GNode<N,E> * n1, GNode<N,E> * n2, Graph<N,E> * g1, Graph<N,E> * g2 )
  {
    checkGraphs(g1, g2);
    if (n1->attr == n2->attr) return .0;
    else  return 3.0 + 2.0 * max_degree;
  }
  
  
  virtual double NodeDeletionCost(GNode<N,E> * n1, Graph<N,E> * g1){ return 1.0; }
  
  virtual double NodeInsertionCost(GNode<N,E> * n2,Graph<N,E> * g2){ return 1.0; }
  
  

  /**
   * Cost of non equal Edge substitution :
   * cost of deleting the egde + cost of re-inserting it + 1
   */
  virtual double EdgeSubstitutionCost(GEdge<int> * e1,GEdge<int> * e2,Graph<N,E> * g1,Graph<N,E> * g2){
    checkGraphs(g1,g2);
    if (e1->attr == e2->attr) return .0;
    else return 3.0;
  }

  
  virtual double EdgeDeletionCost(GEdge<int> * e1,Graph<N,E> * g1){ return 1.0; }
  
  virtual double EdgeInsertionCost(GEdge<int> * e2,Graph<N,E> * g2){ return 1.0; }
};
  


#endif
