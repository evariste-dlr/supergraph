#include <iostream>
#include <fstream>

#include "SymbolicGraph.h"
#include "Dataset.h"

using namespace std;

int main(int argc, char* argv[]){

  if (argc != 4){
    cout << "Usage ! " << argv[0] << "  dortmund_path   N   dataset.json" << endl;
    return 0;
  }
  
  
  string file_A = string(argv[1]) + "A.txt";
  string file_graph_id = string(argv[1]) + "graph_indicator.txt";
  string file_graph_labels = string(argv[1]) + "graph_labels.txt";
  string file_node_labels = string(argv[1]) + "node_labels.txt";
  
  int N = atoi(argv[2]);

  ChemicalDataset<double> dataset;
  ofstream output (argv[2], ofstream::out);


  /* Read graph indicator */

  vector<int> graph_id(N);
  vector<int> node_label(N);
  ifstream ifs_id( file_graph_id, ifstream::in );
  ifstream ifs_nl( file_node_labels, ifstream::in );
  
  int i;
  for (i=0; i<N && ifs_id.good() && ifs_nl.good(); i++){
    ifs_id >> graph_id[i];
    ifs_nl >> node_label[i];
  }
  if (i<N)
    cout << " [W]  EOF before reading all the nodes " << endl;


  /* Create empty graphs */
  int gid = -1;
  int counter = 0;
  graph<int,int>* g = NULL;
  for (int i=0; i<N; i++){
    if (gid != graph_id[i]){
      g = new graph<int,int>();
    }
    
    g->Add(new GNode<int,int>(counter, node_label[i]);
    
    gid = i;
  }


  output << "[" << endl;
