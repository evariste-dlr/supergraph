#include <iostream>
#include <set>

#include "graph.h"
#include "Dataset.h"
#include "SymbolicGraph.h"
#include "ConstantGraphEditDistance.h"
#include "BipartiteGraphEditDistance.h"
#include "MultistartRefinementGraphEditDistance.h"
#include "IPFPGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"
//#include "SupergraphCostFunction.h"
//#include "SuperGraph.h"

#include "include/supergraph.h"
#include "include/utils.h"
#include "include/dataset_json.h"

int main(int argc, char* argv[]){

  std::string ds_path;
  std::string out_graph;
  std::string out_indices;
  std::string out_proj;
  bool project = false;
  float trainprop = 0.1;
  float valprop = 0.1;

  if (argc > 5 && argc < 10){
    ds_path = argv[1];
    trainprop = atof(argv[2]);
    valprop = atof(argv[3]);
    out_graph = argv[4];
    out_indices = argv[5];
    if (trainprop >= 1.0) trainprop = .99;
    else if (trainprop <= .0) trainprop = .01;
    if (valprop >= 1.0 - trainprop) valprop = 1.0 - trainprop - 0.01;
    else if (valprop <= .0) valprop = .01;
  }
  else{
    std::cout << "Usage :  " << argv[0] << "   dataset   trainprop   valprop   outgraph   outindices   [k  [ outproj  p]]] " << std::endl;
    std::cout << "        dataset :  A .ds or .json file listing the graphs" << std::endl;
    std::cout << "      trainprop :  The proportion (in [0;1]) of training examples" << std::endl;
    std::cout << "        valprop :  Proportion (in [0;1]) of validation examples" << std::endl;
    std::cout << "       outgraph :  The outputed supergraph of the train set" << std::endl;
    std::cout << "     outindices :  An output file containing the indices of train, val and test examples" << std::endl;
    std::cout << "        outproj :  Output file with projections on the supergraph " << std::endl;
    std::cout << "              k :  Number of solutions used in the ged estimation" << std::endl;
    std::cout << "              p :  Augmenting data by taking into account at most p% non-optimal projections" << std::endl << std::endl;
    return 0;
  }

  int ksol = 20;     // Nb of solutions for mIPFP
  double paug = .0;  // Percentage of data augmentation
  if (argc > 6)
    ksol = atoi(argv[6]);

  if (argc > 7){
    out_proj = argv[7];
    project = true;
    paug = atof(argv[8]);
  }



  /******************************************
   * Read the dataset
   ******************************************/

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
   * Create the train, validation and test sets
   *****************************************/

  std::default_random_engine randGen;
  randGen.seed(123);
  std::vector<uint> shuf_idx(dataset.size());
  for (uint i=0; i<shuf_idx.size(); i++)  shuf_idx[i] = i;
  std::shuffle( shuf_idx.begin(),  shuf_idx.end(),  randGen );

  std::vector<uint> trainidx(trainprop * dataset.size());
  std::vector<uint> validx(valprop * dataset.size());
  std::vector<uint> testidx(dataset.size() - trainidx.size() - validx.size());

  int N_train = trainidx.size();
  int N_val = validx.size();
  for (uint i=0; i<trainidx.size(); i++)  trainidx[i] = shuf_idx[i];
  for (uint i=0; i<validx.size(); i++)  validx[i] = shuf_idx[N_train+i];
  for (uint i=0; i<testidx.size(); i++)  testidx[i] = shuf_idx[N_train+N_val+i];

  Dataset<int, int, double> trainset; // original dataset for creating the supergraph
  for (uint i=0; i<trainidx.size(); i++){
    trainset.add(new Graph<int,int>(*(dataset[trainidx[i]])), .0);
  }



  /******************************************
   * Initialize the edit distance utilities
   ******************************************/

  struct timeval tv1, tv2;

  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(1000,1,1,2,1,1);
  //CSCostFunction<int,int> * cf = new CSCostFunction<int,int>();
  //EditDistanceCost<std::set<int>, int>* cf = new ChiSquareDistanceCost(1,1,1,1);
  BipartiteGraphEditDistanceMulti<int,int> *ed_init = new BipartiteGraphEditDistanceMulti<int,int>(cf, ksol);
  IPFPGraphEditDistance<int,int> * ipfp = new IPFPGraphEditDistance<int,int>(cf);
  MultistartRefinementGraphEditDistance<int,int> ed(cf, ed_init, ksol, ipfp);



  /******************************************
   * Compute the supergraph in trainset
   ******************************************/

  std::cout << "Computing the supergraph..." << std::endl;

  gettimeofday(&tv1, NULL);
  Graph<int, int>* sg = supergraph<int,int>(trainset, compNode_first, compEdge_first, ed);
  gettimeofday(&tv2, NULL);
  std::cout << std::endl << "--------------" << std::endl;
  std::cout << ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << "s " << std::endl;
  std::cout << "--------------" << std::endl << std::endl;

  std::ofstream ofs_sg (out_graph, std::ofstream::out);
  outputJSON_sg(sg, &ofs_sg);
  ofs_sg.close();


  /******************************************
   * Compute the projections onto the supergraph
   ******************************************/

  if (project){
    //dataset.shuffleize();
    ProjectionDistanceCost<int,int> pcf(ed.getCostFunction(), .001, 1, .0, 1);

    std::cout << "Number of Projections :" << std::endl << std::endl;
    std::list<int*> mappings;
    std::list<Graph<int,int>*> graphs;
    std::list<unsigned int> id_graphs;
    for (int i=0; i<dataset.size(); i++){
      std::list<int*> maps;
      projection(*(dataset[i]), *sg, ed, &maps, &pcf, paug);
      mappings.insert(mappings.end(), maps.begin(), maps.end());
      //std::cout << " graph " << i << " : " << maps.size() << std::endl;
      for (uint k=0; k<maps.size(); k++) {
        graphs.push_back(dataset[i]);
        id_graphs.push_back(i);
      }

      int progress = 51*(float)(i)/dataset.size();
      std::cout << "\r" << "[" << std::string(progress, '=') << std::string(50-progress, ' ') << "]";
      std::cout.flush();
    }
    std::cout << std::endl;

    // update the weight of egdes with projection frequency
    //countFreqProj( sg, mappings, graphs);

    std::ofstream ofs_proj (out_proj, std::ofstream::out);
    outputJSON_proj(sg, mappings, graphs, id_graphs, dataset, &ofs_proj);
    ofs_proj.close();

    for (std::list<int*>::iterator it=mappings.begin(); it!=mappings.end(); it++)
      delete[] *it;
  }


  /******************************************
   * Output the list of train, val and test indices
   ******************************************/

  std::ofstream ofs_idx (out_indices, std::ofstream::out);

  ofs_idx << "{" << std::endl << "  \"train\": [";
  for (uint i=0; i<trainidx.size()-1; i++)
    ofs_idx << trainidx[i] << ", ";

  ofs_idx << trainidx.back() << "]," << std::endl << "  \"val\": [";
  for (uint i=0; i<validx.size()-1; i++)
    ofs_idx << validx[i] << ", ";

  ofs_idx << validx.back() << "]," << std::endl << "  \"test\": [";
  for (uint i=0; i<testidx.size()-1; i++)
    ofs_idx << testidx[i] << ", ";

  ofs_idx << testidx.back() << "]" << std::endl << "}" << std::endl;
  ofs_idx.close();


  delete cf;
  delete ed_init;
  delete ipfp;
  //delete sg; // < deleted by the dataset !

  return 0;
}
