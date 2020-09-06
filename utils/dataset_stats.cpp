#include <iostream>
#include <set>

#include "graph.h"
#include "Dataset.h"
#include "SymbolicGraph.h"

#include "include/dataset_json.h"
#include "include/utils.h"

int main(int argc, char* argv[]){

  std::string ds_path;

  if (argc == 2){
    ds_path = argv[1];
  }
  else{
    std::cout << "Usage :  " << argv[0] << "   dataset   " << std::endl;
    std::cout << "        dataset :  A .ds or .json file listing the graphs" << std::endl << std::endl;
    return 0;
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


  /******************************************
   * Compute stats
   ******************************************/

  int count_nodes = 0;
  int count_edges = 0;
  int positives = 0;
  int negatives = 0;
  for (int i=0; i<dataset.size(); i++){
    count_nodes += dataset[i]->Size();
    count_edges += dataset[i]->getNbEdges();
    if (dataset.getProperty(i) == 0)  negatives ++;
    else positives++;
  }

  const char* _ds_name = strrchr(ds_path.c_str(),'/');
  std::string ds_name;
  if (_ds_name == NULL)  ds_name = ds_path;
  else    ds_name = _ds_name+1;
  std::cout << "\tDataset \t # Graphs \t # nodes \t # Edges \t #+ \t #-" << std::endl;
  std::cout << "\t" << ds_name << " \t " << GREEN << dataset.size() << " \t\t " 
            << float(count_nodes) / dataset.size() << " \t "
            << float(count_edges) / dataset.size() << " \t "
            << positives << " \t " << negatives
            << RESET << std::endl;

  return 0;
}
