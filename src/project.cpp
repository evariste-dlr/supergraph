
/**
 * @file project.cpp
 *
 * project a set of graphs onto a given supergraph.
 */


#include <iostream>
#include <cstdio>
#include <set>

#include "graph.h"
#include "Dataset.h"
#include "SymbolicGraph.h"
#include "BipartiteGraphEditDistance.h"
#include "MultistartRefinementGraphEditDistance.h"
#include "IPFPGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"
//#include "SupergraphCostFunction.h"
//#include "SuperGraph.h"

#include "include/supergraph.h"
#include "include/utils.h"
#include "include/dataset_json.h"

#include "../rapidjson/document.h"
#include "../rapidjson/filereadstream.h"

int main(int argc, char* argv[]){

  std::string ds_path;
  std::string sg_path;
  std::string idx_path;
  std::string out_path;
  if (argc >= 5){
    ds_path = argv[1];
    sg_path = argv[2];
    idx_path = argv[3];
    out_path = argv[4];
  }
  else{
    std::cout << "Usage :  " << argv[0] << "   dataset   supergraph   idx_file   out_file   [k  [p]]  " << std::endl;
    std::cout << "          dataset :  A .ds or .json file listing the graphs to be projected" << std::endl;
    std::cout << "       supergraph :  The prefix of JSON files describing graphs (eg 'sg_' for sg_*.json)" << std::endl;
    std::cout << "         idx_file :  A .json file containing the indices of each subset (one for each supergraph)" << std::endl;
    std::cout << "      output_path :  The output file (JSON)" << std::endl;
    std::cout << "                k :  Number of solutions used in the ged estimation" << std::endl;
    std::cout << "                p :  Data augmentation proportion " << std::endl;
    return 0;
  }

  int ksol = 10;
  if (argc > 5){
    ksol = atoi(argv[5]);
  }

  float paug = 0.3;
  if (argc > 6){
    paug = atof(argv[6]);
  }




  /*****************************************
   * Read the dataset
   *****************************************/
  
  ChemicalDataset<double>* ptr_dataset = NULL;
  const char* ext = strrchr(ds_path.c_str(),'.');
  if (ext != NULL && strcmp(ext, ".json") == 0){
    ptr_dataset = new ChemicalDataset<double>();
    dataset_json(argv[1], *ptr_dataset);
  }
  else
    ptr_dataset = new ChemicalDataset<double>(argv[1]);
  
  ChemicalDataset<double>& dataset = *ptr_dataset;
  dataset.shuffleize();
    


  /*****************************************
   * Read the indices file
   *****************************************/
  
  FILE * fp = fopen(idx_path.c_str(), "rb");
  char readBuffer[65536];
  rapidjson::FileReadStream ifs(fp, readBuffer, sizeof(readBuffer));
  rapidjson::Document d;
  d.ParseStream(ifs);
  if (d.HasParseError()){
    std::cout << "[E]  Parse error in the supergraph JSON file" << std::endl;
    fclose(fp);
    return 1;
  }

  
  /******************************************
   * Initialize the edit distance utilities
   ******************************************/
  
  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(1000,1,1,3,1,1);
  //CSCostFunction<int,int> * cf = new CSCostFunction<int,int>();
  //EditDistanceCost<std::set<int>, int>* cf = new ChiSquareDistanceCost(1,1,1,1);
  BipartiteGraphEditDistanceMulti<int,int> *ed_init = new BipartiteGraphEditDistanceMulti<int,int>(cf, ksol);
  IPFPGraphEditDistance<int,int> * ipfp = new IPFPGraphEditDistance<int,int>(cf);
  ipfp->compute_equiv_mappings(true);
  MultistartRefinementGraphEditDistance<int,int> ed(cf, ed_init, ksol, ipfp);
    
    
  

  /***** Number of supergraphs *****/
  int N_SG = 10;


  /*** For Each Supergraph ***/
  for (int i=0; i<N_SG; i++){


  /*****************************************
   * Read the supergraph
   *****************************************/

    std::string file_path = sg_path+itoa(i,10)+".json";
    FILE * fp = fopen(sg_path.c_str(), "rb");
    //char readBuffer[65536];
    rapidjson::FileReadStream ifs(fp, readBuffer, sizeof(readBuffer));
    rapidjson::Document d;
    d.ParseStream(ifs);
    if (d.HasParseError()){
      std::cout << "[E]  Parse error in the supergraph JSON file" << std::endl;
      fclose(fp);
      return 1;
    }
    
    // Non-directed supergraph
    Graph<int, int>* sg = new Graph<int,int>(false);
    
    assert(d.IsObject());
    assert(d.HasMember("nodes"));
    assert(d.HasMember("links"));
    
    const rapidjson::Value& nodes = d["nodes"];
    assert(nodes.IsArray());
    for (rapidjson::Value::ConstValueIterator itr = nodes.Begin(); itr != nodes.End(); ++itr){
      int id = (*itr)["id"].GetInt();
      int atom = (*itr)["atom"].GetInt();
      GNode<int,int>* n = new GNode<int,int>(id, atom);
      sg->Add(n);
    }
    
    const rapidjson::Value& edges = d["links"];
    assert(edges.IsArray());
    for (rapidjson::Value::ConstValueIterator itr = edges.Begin(); itr != edges.End(); ++itr){
      int from = (*itr)["source"].GetInt();
      int to   = (*itr)["target"].GetInt();
      sg->Link(from, to, 1);
    }
    
    fclose(fp);


  

  }
  
  std::ofstream output (out_path, std::ofstream::out);
  outputJSON_proj(sg, mappings, graphs, id_graphs, cdataset, &output);
  output.close();
  

  for (std::list<int*>::iterator it=mappings.begin(); it!=mappings.end(); it++)
    delete *it;
  
  delete cf;
  delete ed_init;
  delete ipfp;
  delete sg;
  
  return 0;
}
