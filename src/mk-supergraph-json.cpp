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
#include "SupergraphCostFunction.h"
#include "SuperGraph.h"

#include "supergraph.h"
#include "utils.h"

int main(int argc, char* argv[]){

  std::string train_path;
  std::string test_path;
  std::string out_path;
  std::string val_path;
  bool validation=false;
  if (argc == 4){
    train_path = argv[1];
    test_path = argv[2];
    out_path = argv[3];
  }
  else if (argc > 4){
    train_path = argv[1];
    test_path = argv[2];
    val_path  = argv[3];
    out_path = argv[4];
    validation=true;
  }
  else{
    std::cout << "Usage :  " << argv[0] << "   trainset  testset valset output_path  [k  [p]] " << std::endl;
    std::cout << "         trainset :  A .ds file listing the training graph files" << std::endl;
    std::cout << "          testset :  A .ds file listing the test examples" << std::endl;
    std::cout << "           valset :  A .ds file listing the validation examples" << std::endl;
    std::cout << "      output_path :  The output directory" << std::endl;
    std::cout << "                k :  Number of solutions used in the ged estimation (default 10)" << std::endl;
    std::cout << "                p :  Augmenting data by taking into account at most p% non-optimal projections (default 0)" << std::endl;
    return 0;
  }

  int ksol = 10;     // Nb of solutions for mIPFP
  double paug = 0.0; // Percentage of data augmentation
  if (argc > 5)
    ksol = atoi(argv[5]);

  if (argc > 6)
    paug = atof(argv[6]);


  struct timeval tv1, tv2;
  
  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(1000,1,1,2,1,1);
  //CSCostFunction<int,int> * cf = new CSCostFunction<int,int>();
  //EditDistanceCost<std::set<int>, int>* cf = new ChiSquareDistanceCost(1,1,1,1);
  BipartiteGraphEditDistanceMulti<int,int> *ed_init = new BipartiteGraphEditDistanceMulti<int,int>(cf, ksol);
  IPFPGraphEditDistance<int,int> * ipfp = new IPFPGraphEditDistance<int,int>(cf);
  MultistartRefinementGraphEditDistance<int,int> ed(cf, ed_init, ksol, ipfp);
  
  ChemicalDataset<double> ctrainset(train_path.c_str()); // original dataset for creating the supergraph
  ChemicalDataset<double> ctestset(test_path.c_str());   // Dataset used for computing tests projections
  Dataset<int, int, double> dataset;
  for (int i=0; i<ctrainset.size(); i++){
    dataset.add(new Graph<int,int>(*(ctrainset[i])), .0);
  }

  gettimeofday(&tv1, NULL);
  Graph<int, int>* sg = supergraph<int,int>(dataset, compNode_first, compEdge_first, ed);
  gettimeofday(&tv2, NULL);
  std::cout << std::endl << "--------------" << std::endl;
  std::cout << ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << "s " << std::endl;
  std::cout << "--------------" << std::endl;
  
  
  // Tester les projections de chaque graphe de la base sur le supergraphe
  ctrainset.shuffleize();
  ProjectionDistanceCost<int,int> pcf(ed.getCostFunction(), .001, 1, .0, 1);
  //std::cout << "id \t size \t dist \t success" << std::endl << "--------------------------------" << std::endl;
  std::cout << "Train projections" << std::endl;
  std::list<int*> mappings;
  std::list<Graph<int,int>*> graphs;
  std::list<unsigned int> id_graphs;
  for (int i=0; i<ctrainset.size(); i++){
    //int* mapping = new int[ctrainset[i]->Size()];
    std::list<int*> maps;
    projection(*(ctrainset[i]), *sg, ed, &maps, &pcf, paug);
    mappings.insert(mappings.end(), maps.begin(), maps.end());
    std::cout << "mappings : " << maps.size() << std::endl;
    for (uint k=0; k<maps.size(); k++) {
      graphs.push_back(ctrainset[i]);
      id_graphs.push_back(i);
    }
  }

  std::cout << "Test projections" << std::endl;
  std::list<int*> map_test;
  std::list<Graph<int,int>*> graphs_test;
  std::list<unsigned int> id_graphs_test;
  for (int i=0; i<ctestset.size(); i++){
    //int* mapping = new int[ctrainset[i]->Size()];
    std::list<int*> maps;
    projection(*(ctestset[i]), *sg, ed, &maps, &pcf, paug);
    map_test.insert(map_test.end(), maps.begin(), maps.end());
    std::cout << "mappings : " << maps.size() << std::endl;
    for (uint k=0; k<maps.size(); k++) {
      graphs_test.push_back(ctestset[i]);
      id_graphs_test.push_back(i);
    }
  }

  // If there is a validation set
  if (validation){
    ChemicalDataset<double> cvalset(val_path.c_str());   // Dataset used for computing validation projections
    std::cout << "Valid. projections" << std::endl;
    std::list<int*> map_val;
    std::list<Graph<int,int>*> graphs_val;
    std::list<unsigned int> id_graphs_val;
    for (int i=0; i<cvalset.size(); i++){
      //int* mapping = new int[ctrainset[i]->Size()];
      std::list<int*> maps;
      projection(*(cvalset[i]), *sg, ed, &maps, &pcf, paug);
      map_val.insert(map_val.end(), maps.begin(), maps.end());
      std::cout << "mappings : " << maps.size() << std::endl;
      for (uint k=0; k<maps.size(); k++) {
	graphs_val.push_back(cvalset[i]);
	id_graphs_val.push_back(i);
      }
    }
    std::ofstream out_val (out_path+"/val.json", std::ofstream::out); // val projections
    outputJSON_proj(sg, map_val, graphs_val, id_graphs_val, cvalset, &out_val);
    out_val.close();
  }

  // update the weight of egdes with projection frequency
  //countFreqProj( sg, mappings, graphs);

  std::ofstream out_sg (out_path+"/supergraph.json", std::ofstream::out);
  std::ofstream out_tr (out_path+"/train.json", std::ofstream::out); // train projections
  std::ofstream out_ts (out_path+"/test.json", std::ofstream::out); // test projections
  
  
  outputJSON_sg(sg, &out_sg);
  outputJSON_proj(sg, mappings, graphs, id_graphs, ctrainset, &out_tr);
  outputJSON_proj(sg, map_test, graphs_test, id_graphs_test, ctestset, &out_ts);
  out_sg.close();
  out_tr.close();
  out_ts.close();
  

  for (std::list<int*>::iterator it=mappings.begin(); it!=mappings.end(); it++)
    delete[] *it;
  
  delete cf;
  delete ed_init;
  delete ipfp;
  //delete sg; // < deleted by the dataset !
  
  return 0;
}
