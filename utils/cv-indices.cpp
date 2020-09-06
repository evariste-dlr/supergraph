
/**
 * @file cv-indices.cpp
 *
 *  Create indices files for cross validation from a dataset.
 *  The program outputs a json file containing the train, valid and test indices
 */


#include <iostream>
//#include <set>

#include "graph.h"
#include "Dataset.h"
#include "SymbolicGraph.h"

#include "include/utils.h"
#include "include/dataset_json.h"

int main(int argc, char* argv[]){

  std::string ds_path;
  std::string out_indices;
  int k_fold = 10;
  float trainprop = 0.3;
  float valprop = 0.3;

  if (argc == 5){
    ds_path = argv[1];
    k_fold = atoi(argv[2]);
    trainprop = atof(argv[3]);
    out_indices = argv[4];
    if (trainprop >= 1.0) trainprop = .99;
    else if (trainprop <= .0) trainprop = .01;
    if (valprop >= 1.0 - trainprop) valprop = 1.0 - trainprop - 0.01;
    else if (valprop <= .0) valprop = .01;
  }
  else{
    std::cout << "Usage :  " << argv[0] << "   dataset   nb-cross   trainprop   out_name " << std::endl;
    std::cout << "        dataset :  A .ds or .json file listing the graphs" << std::endl;
    std::cout << "       nb-cross :  Number k of cross-validation test sets" << std::endl;
    std::cout << "      trainprop :  The proportion (in [0;1]) of training examples in the remaining set" << std::endl;
    std::cout << "       out_name :  The outputed json files : 1, named out_name_k.json, per different testset k" << std::endl;
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




  /*****************************************
   * Shuffle the dataset and create k subsets for cross-validation
   *****************************************/

  std::default_random_engine randGen;
  randGen.seed(123);
  std::vector<uint> shuf_idx(dataset.size());
  for (uint i=0; i<shuf_idx.size(); i++)  shuf_idx[i] = i;
  std::shuffle( shuf_idx.begin(),  shuf_idx.end(),  randGen );

  std::vector<uint> testidx(dataset.size() / k_fold);
  std::vector<uint> trainidx(trainprop * (dataset.size() - testidx.size()));
  std::vector<uint> validx(dataset.size() - trainidx.size() - testidx.size());

  int N_train = trainidx.size();
  int N_test = testidx.size();


  /*      For each test set      */
  for (int k_test=0; k_test < k_fold;  k_test++){

  for (int i=0; i<N_test; i++)  testidx[i] = shuf_idx[k_test * N_test + i];
  int j = k_test == 0 ? N_test : 0;
  for (int i=0; j<dataset.size(); i++,(j==k_test*N_test-1) ? j+=N_test+1 : j++){
    if (i < N_train)
      trainidx[i] = shuf_idx[j];
    else if (i-N_train < validx.size())
      validx[i-N_train] = shuf_idx[j];
  }

  Dataset<int, int, double> trainset; // original dataset for creating the supergraph
  for (uint i=0; i<trainidx.size(); i++){
    trainset.add(new Graph<int,int>(*(dataset[trainidx[i]])), .0);
  }



  /******************************************
   * Output the list of train, val and test indices
   ******************************************/

  std::ofstream ofs_idx (out_indices+"_"+std::to_string(k_test)+".json", std::ofstream::out);

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


  }

  delete ptr_dataset;

  return 0;
}
