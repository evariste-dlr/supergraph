#include <iostream>
#include <iomanip>
#include <bitset>

#include "include/stars.h"
#include "include/utils.h"
#include "include/dataset_json.h"

// graph-lib includes
#include "SymbolicGraph.h"
#include "ConstantGraphEditDistance.h"
#include "BipartiteGraphEditDistance.h"
#include "MultistartRefinementGraphEditDistance.h"
#include "IPFPGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"

using namespace std;


int main(int argc, char* argv[]){

  if (argc < 3){
    cout << "Usage : " << argv[0] << "  dataset_file  output_file  [ds_json]" << endl;
    cout << "      dataset_file   .ds file or .json file" << endl;
    cout << "      output_file    .json file" << endl;
    cout << "      ds_json=0      1 for a json dataset file" << endl;
    return 0;
  }

  bool ds_json = false;
  if (argc == 4 && atoi(argv[3]) == 1)
    ds_json = true;

  cout << GREEN << "[I] Loading dataset..." << RESET << endl;
  ChemicalDataset<double>* ptr_dataset = NULL;
  if (ds_json){
    ptr_dataset = new ChemicalDataset<double>();
    dataset_json(argv[1], *ptr_dataset);
  }
  else
    ptr_dataset = new ChemicalDataset<double>(argv[1]);

  ChemicalDataset<double>& dataset = *ptr_dataset;

  cout << GREEN << "[I] Extracting all motifs of the dataset..." << RESET << endl;

  vector<unordered_map<uint,uint> > * distrib = new vector<unordered_map<uint,uint> >[dataset.size()];
  for (int i=0; i<dataset.size(); i++)
    distrib[i].resize(dataset[i]->Size());

  unordered_map<std::string, uint>* dict = star_dictionnary(dataset, distrib);


  // Résultats
  ofstream output (argv[2], ofstream::out);
  output << "[" << endl;
  cout << endl;
  cout << "---------------------------------------------------" << endl;
  cout << "  Dictionnary " << endl;
  cout << "                   code (40 first chars)        idx" << endl;
  cout << "---------------------------------------------------" << endl;
  for (unordered_map<std::string,uint>::iterator it=dict->begin(); it!=dict->end(); it++){
    cout << std::setw(40) << it->first;
    cout << " \t " << it->second << endl;
  }

  for (int idx=0; idx<dataset.size(); idx++){
          //cout << endl << "--------------------" << endl << "   Graph " << idx << endl;
          //cout << "Node \t motifs" << endl;
          output << "  {" << endl;
          output << "    \"graph\": " << idx << "," << endl;
          output << "    \"motifs\": [" << endl;
          for (int i=0; i<dataset[idx]->Size(); i++){
            //cout << i << "\t ";
            output << "      [";
            for (uint m=0; m<dict->size(); m++){
              if (distrib[idx][i].count(m) != 0){
                output << distrib[idx][i][m];
              }
              else{
                output << "0";
              }
              if (m < dict->size()-1) output << ",";
            }
            if (i < dataset[idx]->Size()-1) output << "]," << endl;
            else  output << "]" << endl;
          }
          output << "    ]" << endl;
          if (idx < dataset.size()-1)  output << "  }," << endl;
          else  output << "  }" << endl;
  }

  output << "]" << endl;
  output.close();

  delete[] distrib;
  delete dict;
  delete ptr_dataset;

  return 0;
}
