/**
 * @file utils.h  Utilities
 */

#ifndef __GCNN_UTILS_H__
#define __GCNN_UTILS_H__

#include <iostream>
#include <string>
#include <set>
#include <list>

#include "graph.h"
#include "ECMapping.h"



const std::string RESET         = "\x1b[0m";
const std::string BLINK         = "\x1b[5m";
const std::string GREEN         = "\x1b[92m";
const std::string GREEN_BACK    = "\x1b[42m\x1b[1m";
const std::string RED           = "\x1b[91m\x1b[1m";
const std::string RED_BACK      = "\x1b[101m";
const std::string CYAN          = "\x1b[96m";


/**
 * A set as a string
 */
std::string setToString(const std::set<int>& set ){
  std::string str("");
  for (std::set<int>::iterator it = set.begin(); it!=set.end(); it++){
    str +=  std::to_string(*it) + ";";
  }
  return str;
}


void outputGexf_set( Graph<std::set<int>,int>* sg,
		     std::list<int*>& mappings,
		     Dataset<int,int,double>& dataset,
		     std::ostream * out = &(std::cout)){
  // GEXF Format
  *out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  *out << "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">" << std::endl;
  *out << "  <graph mode=\"static\" defaultedgetype=\"undirected\">" << std::endl;
  *out << "    <attributes class=\"node\"> <attribute id=\"0\" title=\"atom\" type=\"string\"/>" << std::endl;
  for (unsigned int i=1; i<=mappings.size(); i++){
    *out << "      <attribute id=\""<< i << "\" title=\"proj_" << i << "\" type=\"integer\"/>" << std::endl;
  }
  *out << "    </attributes>" << std::endl;
  *out << "    <nodes> " << std::endl;
  for (int i=0; i<sg->Size(); i++){
    *out << "     <node id=\"" << i+1 << "\">" << std::endl;
    *out << "       <attvalues><attvalue for=\"0\" value=\"" << setToString((*sg)[i]->attr) << "\" />" << std::endl;
    int g=0;
    for (std::list<int*>::iterator it=mappings.begin(); it!=mappings.end(); it++){
      int k=0;
      while(k<dataset[g]->Size() && (*it)[k]!=i) k++;
      if (k < dataset[g]->Size())
	*out << "       <attvalue for=\"" << g+1 << "\" value=\""<<1<<"\" />" << std::endl;
      else
	*out << "       <attvalue for=\"" << g+1 << "\" value=\""<<0<<"\" />" << std::endl;
      g++;
    }
    *out << "       </attvalues>" << std::endl << "     </node>" << std::endl;
  }
  *out << "    </nodes> <edges>" << std::endl;

  for (int i=0; i<sg->Size()-1; i++){
    for (int j=i+1; j<sg->Size(); j++){
      if (sg->isLinked(i,j)){
	*out << "     <edge source=\"" << i+1 << "\" target=\"" << j+1 << "\" />" << std::endl;
	//*out << "    label " << sg->getEdge(i,j)->attr << std::endl;
      }
    }
  }
  *out << "</edges></graph> </gexf>" << std::endl;

}


void outputGexf_int( Graph<int,int>* sg,
		     std::list<int*>& mappings,
		     std::list<Graph<int,int>* >&  graphs,
		     Dataset<int, int, double>& dataset,
		     std::ostream * out = &(std::cout) )
{
  // GEXF Format
  *out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  *out << "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">" << std::endl;
  *out << "  <graph mode=\"static\" defaultedgetype=\"undirected\">" << std::endl;
  *out << "    <attributes class=\"node\"> "<< std::endl;
  *out << "      <attribute id=\"0\" title=\"atom\" type=\"string\"/>" << std::endl;
  for (unsigned int i=1; i<=mappings.size(); i++){
    *out << "      <attribute id=\""<< i << "\" title=\"proj_" << i << "\" type=\"integer\"/>" << std::endl;
  }
  *out << "    </attributes>" << std::endl;
  *out << "    <attributes class=\"edge\"> "<< std::endl;
  for (unsigned int i=1; i<=mappings.size(); i++){
    *out << "      <attribute id=\""<< i << "\" title=\"proj_" << i << "\" type=\"integer\"/>" << std::endl;
  }
  *out << "    </attributes>" << std::endl;
  *out << "    <nodes> " << std::endl;
  for (int i=0; i<sg->Size(); i++){
    *out << "     <node id=\"" << i+1 << "\">" << std::endl;
    *out << "       <attvalues><attvalue for=\"0\" value=\"" << (*sg)[i]->attr << "\" />" << std::endl;
    int g=0;
    std::list<Graph<int,int>*>::iterator it_graphs = graphs.begin();
    for (std::list<int*>::iterator it=mappings.begin(); it!=mappings.end(); it++){
      int k=0;
      while(k<(*it_graphs)->Size() && (*it)[k]!=i) k++;
      if (k < (*it_graphs)->Size())
	*out << "       <attvalue for=\"" << g+1 << "\" value=\""<<1<<"\" />" << std::endl;
      else
	*out << "       <attvalue for=\"" << g+1 << "\" value=\""<<0<<"\" />" << std::endl;
      it_graphs++; g++;
    }
    *out << "       </attvalues>" << std::endl << "     </node>" << std::endl;
  }
  *out << "    </nodes> <edges>" << std::endl;

  for (int i=0; i<sg->Size()-1; i++){
    for (int j=i+1; j<sg->Size(); j++){
      if (sg->isLinked(i,j)){
	*out << "     <edge source=\"" << i+1 << "\" target=\"" << j+1 << "\">" << std::endl;
	int g=0;
	std::list<Graph<int,int>*>::iterator it_graphs = graphs.begin();
	for (std::list<int*>::iterator it=mappings.begin(); it!=mappings.end(); it++){
	  int k=0;  while(k<(*it_graphs)->Size() && (*it)[k]!=i) k++;
	  int ri = k;
	  k=0;  while(k<(*it_graphs)->Size() && (*it)[k]!=j) k++;
	  int rj = k;
	  if (ri < (*it_graphs)->Size() && rj < (*it_graphs)->Size() && (*it_graphs)->isLinked(ri,rj)){
	    *out << "       <attvalue for=\"" << g+1 << "\" value=\""<<1<<"\" />" << std::endl;
	  }
	  else
	    *out << "       <attvalue for=\"" << g+1 << "\" value=\""<<0<<"\" />" << std::endl;
	  g++;
	  it_graphs++;
	}
	*out << "     </edge>" << std::endl;
	//*out << "    label " << sg->getEdge(i,j)->attr << std::endl;
      }
    }
  }
  *out << "</edges></graph> </gexf>" << std::endl;

}



/**
 *  Output a supergraph in JSON format
 */
void outputJSON_sg( Graph<int,int>* sg,
		    std::ostream * out_sg = &(std::cout))
{
  // Supergraph nodes
  *out_sg << "{" << std::endl << " \"nodes\" : [" << std::endl;
  for (int i=0; i<sg->Size(); i++){
    *out_sg << "  {\"id\": " << i << ", \"atom\": " << (*sg)[i]->attr;
    if (i < sg->Size()-1) *out_sg << "}," << std::endl;
    else *out_sg << "}" << std::endl;
  }
  *out_sg << " ]," << std::endl << "\"links\" : [" << std::endl;

  // Supergraph edges
  bool first_edge = true;
  for (int i=0; i<sg->Size()-1; i++){
    for (int j=i+1; j<sg->Size(); j++){
      if (sg->isLinked(i,j)){
	if (first_edge)  first_edge = false;
	else *out_sg << "," << std::endl;
	*out_sg << "  {\"source\": " << i
	        << ", \"target\": " << j
	        << ", \"bound\": " << sg->getEdge(i,j)->attr
	        << "}";
      }
    }
  }

  *out_sg << std::endl << " ]" << std::endl << "}" << std::endl;
}


/**
 * Output projections on a supergraph in JSON format
 */
void outputJSON_proj( Graph<int,int>* sg,
		     std::list<ECMapping*>& mappings,
		     std::list<Graph<int,int>* >&  graphs,
		     std::list<unsigned int>& id_graphs,
		     Dataset<int, int, double>& dataset,
		     std::ostream * out_proj = &(std::cout))
{
  bool first_edge = true;

  // Projections
  int n_graph = 0;
  *out_proj << "[" << std::endl;
  std::list<Graph<int,int>*>::iterator it_graphs = graphs.begin();
  std::list<unsigned int>::iterator it_id = id_graphs.begin();
  for (std::list<ECMapping*>::iterator it=mappings.begin(); it!=mappings.end(); it++){
    // Nodes
    *out_proj << " {" << std::endl << "  \"nodes\": [";
    for (int k=0; k<(*it_graphs)->Size(); k++){
      *out_proj << (**it)[k];
      if (k<(*it_graphs)->Size()-1) *out_proj << ", ";
    }

    // Edges
    *out_proj << "]," << std::endl << "  \"links\": [";
    first_edge=true;
    for (int k=0; k<(*it_graphs)->Size(); k++){
      for (int l=k+1; l<(*it_graphs)->Size(); l++){
        if ((**it)[k] < sg->Size() && (**it)[l] < sg->Size()){
	 if (sg->isLinked((**it)[k], (**it)[l]) && (*it_graphs)->isLinked(k,l)){
	  if (first_edge){
	    *out_proj << "[";
	    first_edge = false;
	  }else *out_proj << ", [";
	  *out_proj  << k << "," << l << "]";
	 }
	}
      }
    }
    *out_proj << "]," << std::endl;

    // Graph id
    if (*(std::prev(it_graphs)) != *it_graphs) n_graph++;
    *out_proj << "  \"graph\": " << n_graph << "," << std::endl;

    // Value (class or regression value)
    *out_proj << "  \"value\": " << dataset.getProperty(*it_id) << std::endl << " }";
    if (std::next(it) != mappings.end()) *out_proj << ",";
    *out_proj << std::endl;


    it_graphs++; it_id++;
  }
  *out_proj << "]" << std::endl;

}



/**
 * @brief Count the frequencies of each edge of a graph in the given projections
 */
void countFreqProj( Graph<int,int> * sg,
		    std::list<int*>& proj,
		    std::list<Graph<int,int>* >& graphs )
{
  // For each edge e in the graph
  for (int i=0; i<sg->Size(); i++){
    for (int j=0; j<sg->Size(); j++){
      GEdge<int> * e = sg->getEdge(i,j);
      if (e != NULL){
	// count the total number of projections containing e
	int count = 0;
	std::list<Graph<int,int>*>::iterator it_graphs = graphs.begin();
	for (std::list<int*>::iterator it=proj.begin(); it!=proj.end(); it++){
	  for (int k=0; k<(*it_graphs)->Size(); k++){
	    if ((*it)[k] == i){
	      for (int l=0; l<(*it_graphs)->Size(); l++){
		if (k != l &&  (*it)[l] == j && (*it_graphs)->isLinked(k,l)){
		  count++;
		  break;
		}
	      }
	    }
	  }
	  it_graphs++;
	}
	e->attr = count;
      }
    }
  }

}

#endif // __UTILS_H__
