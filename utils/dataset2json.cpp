#include <iostream>
#include "Dataset.h"

using namespace std;

int main(int argc, char* argv[]){

  if (argc != 3){
    cout << "Usage ! " << argv[0] << "  dataset.ds  dataset.json" << endl;
    return 0;
  }

  ChemicalDataset<double> dataset(argv[1]);
  ofstream output (argv[2], ofstream::out);

  output << "[" << endl;

  for (int i=0; i<dataset.size(); i++){
    Graph<int,int> * g = dataset[i];
    output << "  {" << endl;
    output << "    \"graph\": " << i << "," << endl;
    output << "    \"value\": " << dataset.getProperty(i) << "," << endl;

    output << "    \"nodes\": [";
    for (int v=0; v<g->Size(); v++){
      output << (*g)[v]->attr;
      if (v<g->Size()-1) output << ",";
    }
    output << "]," << endl;

    output << "    \"adj\": [" << endl;
    for (int v=0; v<g->Size(); v++){
      output << "      [";
      for (int w=0; w<g->Size(); w++){
        if (g->isLinked(v,w))
          output << g->getEdge(v,w)->attr;
        else
          output << "0";
        if (w < g->Size()-1) output << ",";
      }
      if (v < g->Size()-1)  output << "]," << endl;
      else  output << "]" << endl;
    }
    output << "    ]" << endl;
    if (i < dataset.size()-1) output << "  }," << endl;
    else  output << "  }" << endl;
  }

  output << "]" << endl;
  output.close();

  return 0;
}
