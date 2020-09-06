
/**
 * @file compute-supergraph-cv.cpp
 *
 *  k-fold cross-validation version of compute-supergraph.
 *  This version produces k * 3 files, one set for each
 *  selected testset :
 *  * supergraph_k.json  :  A supergraph of the training set
 *  * indices_k.json  :     The indices in the dataset for training, validation and test sets
 *  * project_k.json  :     The projections of each graph onto the supergraph
 */


#include <iostream>
#include <iomanip>
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

#define VERBOSE

#include "include/supergraph.h"
#include "include/utils.h"
#include "include/dataset_json.h"

int main(int argc, char* argv[]){

  std::string ds_path;
  std::string outdir;
  int k_fold = 1;
  float trainprop = 0.3;
  float valprop = 0.3;
  bool use_idFile = false;
  std::string idx_file;

  if (argc > 3 && argc < 8){
    ds_path = argv[1];
    const char* ext_ = strrchr(argv[2],'.');
    if (ext_ != NULL && strcmp(ext_, ".json") == 0){
      idx_file = argv[2];
      use_idFile = true;
      outdir = argv[3];
    }
    else{
      k_fold = atoi(argv[2]);
      trainprop = atof(argv[3]);
      outdir = argv[4];
      if (trainprop >= 1.0) trainprop = .99;
      else if (trainprop <= .0) trainprop = .01;
      if (valprop >= 1.0 - trainprop) valprop = 1.0 - trainprop - 0.01;
      else if (valprop <= .0) valprop = .01;
    }
  }
  else{
    std::cout << "Usage :  " << argv[0] << "   dataset   {idx_file  |  nb-cross  trainprop}   outdir   [k  [p]] " << std::endl;
    std::cout << "        dataset :  A .ds or .json file listing the graphs" << std::endl;
    std::cout << "       idx_file :  A .json file containing the train, test and valid list of indices" << std::endl;
    std::cout << "       nb-cross :  Number k of cross-validation test sets" << std::endl;
    std::cout << "      trainprop :  The proportion (in [0;1]) of training examples in the remaining set" << std::endl;
    std::cout << "         outdir :  Output directory, k * 3 files will be outputed" << std::endl;
    std::cout << "              k :  Number of solutions used in the ged estimation" << std::endl;
    std::cout << "              p :  Augmenting data by taking into account at most p% non-optimal projections" << std::endl;
    std::cout << "                   Negative value <=> do not compute projections" << std::endl << std::endl;
    return 0;
  }

  int ksol = 20;     // Nb of solutions for mIPFP
  double paug = .0;  // Percentage of data augmentation
  if (use_idFile){
    if (argc > 4)
      ksol = atoi(argv[4]);
    if (argc > 5)
      paug = atof(argv[5]);
  }
  else{
    if (argc > 5)
      ksol = atoi(argv[5]);
    if (argc > 6)
      paug = atof(argv[6]);
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
   * Initialize the edit distance utilities
   ******************************************/

  struct timeval tv1, tv2;

  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(1000,1,1,3,1,1);
  //CSCostFunction<int,int> * cf = new CSCostFunction<int,int>();
  //EditDistanceCost<std::set<int>, int>* cf = new ChiSquareDistanceCost(1,1,1,1);
  BipartiteGraphEditDistanceMulti<int,int> *ed_init = new BipartiteGraphEditDistanceMulti<int,int>(cf, ksol);
  IPFPGraphEditDistance<int,int> * ipfp = new IPFPGraphEditDistance<int,int>(cf);
  ipfp->compute_equiv_mappings(false);
  MultistartRefinementGraphEditDistance<int,int> ed(cf, ed_init, ksol, ipfp);

  int N_train = 0;
  int N_test = 0;

  if (!use_idFile) {
    N_test = dataset.size() / k_fold;
    N_train = trainprop * (dataset.size() - N_test);
  }

    /*****************************************
     * Shuffle the dataset and create k subsets for cross-validation
     *****************************************/

  std::default_random_engine randGen;
  randGen.seed(123);
  std::vector<uint> shuf_idx(dataset.size());
  for (uint i=0; i<shuf_idx.size(); i++)  shuf_idx[i] = i;
  std::shuffle( shuf_idx.begin(),  shuf_idx.end(),  randGen );

  std::vector<uint> testidx(N_test);
  std::vector<uint> trainidx(N_train);
  std::vector<uint> validx(dataset.size() - trainidx.size() - testidx.size());


  /*      For each test set      */
  //for (int k_test=0; k_test < k_fold;  k_test++){
  for (int k_test=0; k_test < 1;  k_test++){

    if (use_idFile){
      load_subsets(idx_file, trainidx, testidx, validx);
    }
    else{
      for (int i=0; i<N_test; i++)  testidx[i] = shuf_idx[k_test * N_test + i];
      int j = k_test == 0 ? N_test : 0;
      for (int i=0; j<dataset.size(); i++,(j==k_test*N_test-1) ? j+=N_test+1 : j++){
        if (i < N_train)
          trainidx[i] = shuf_idx[j];
        else if (i-N_train < validx.size())
          validx[i-N_train] = shuf_idx[j];
      }
    }

    Dataset<int, int, double> trainset; // original dataset for creating the supergraph
    for (uint i=0; i<trainidx.size(); i++){
      trainset.add(new Graph<int,int>(*(dataset[trainidx[i]])), .0);
    }


    std::string out_graph(outdir+"/sg_");
    std::string out_proj(outdir+"/proj_");
    std::string out_indices(outdir+"/idx_");
    out_graph += std::to_string(k_test) + ".json";
    out_proj += std::to_string(k_test) + ".json";
    out_indices += std::to_string(k_test) + ".json";



    /******************************************
     * Compute the supergraph of the trainset
     ******************************************/

    std::cout << std::endl << std::endl;
    std::cout << CYAN << "  [Supergraph]  " << RESET << "(" << k_test+1 << "/"<<k_fold <<")" << std::endl;

    #ifdef BIPARTITE_MATRIX
      bool use_bipartite = true;
    #else
      bool use_bipartite = false;
    #endif

    gettimeofday(&tv1, NULL);
    Graph<int, int>* sg = supergraph<int,int>(trainset, compNode_first, compEdge_first, ed, ksol, use_bipartite);
    gettimeofday(&tv2, NULL);

    std::cout << "  [" << ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << "s " << "]" << std::endl;


    std::ofstream ofs_sg (out_graph, std::ofstream::out);
    outputJSON_sg(sg, &ofs_sg);
    ofs_sg.close();




    /******************************************
     * Compute the projections onto the supergraph
     ******************************************/
    
    if (paug >= 0){

      // Activate the computation of all equivalent mappings
      ipfp->compute_equiv_mappings(true);

      std::cout << CYAN << "  [Projections]" << RESET << std::endl;
      std::cout << std::fixed;
      std::cout << std::setprecision(0);

      std::list<ECMapping*> mappings;
      std::list<Graph<int,int>*> graphs;
      std::list<unsigned int> id_graphs;

      #if NESTED_TYPE != 0
        int num_group_threads = NUM_THREADS/ksol;
        #ifdef BIPARTITE_MATRIX
          num_group_threads = NUM_THREADS;
        #endif
        std::vector<MultistartRefinementGraphEditDistance<int,int> *> local_ged(num_group_threads);
        for (unsigned int _g=0; _g<local_ged.size(); _g++)
          local_ged[_g] = ed.clone();

        #pragma omp parallel num_threads(num_group_threads)
        {
      #endif

      ProjectionDistanceCost<int,int> pcf(ed.getCostFunction(), .001, 1, .0, 1);


      #if NESTED_TYPE != 0
        #pragma omp for ordered schedule(dynamic,1)
      #endif
      for (int i=0; i<dataset.size(); i++){

        std::list<ECMapping*> maps;
        #if NESTED_TYPE != 0
          MultistartRefinementGraphEditDistance<int,int> * lged = local_ged[omp_get_thread_num()];
        #else
          MultistartRefinementGraphEditDistance<int,int> * lged = &ed;
        #endif

        projection(*(dataset[i]), *sg, *lged, maps, &pcf, paug);

        #if NESTED_TYPE != 0
          #pragma omp ordered
          {
        #endif
	mappings.insert(mappings.end(), maps.begin(), maps.end());
        //std::cout << " graph " << i << " : " << maps.size() << std::endl;
        for (uint k=0; k<maps.size(); k++) {
          graphs.push_back(dataset[i]);
          id_graphs.push_back(i);
        }
        #if NESTED_TYPE != 0
      } // end omp ordered
          if(omp_get_thread_num() == 0){
        #endif

        float prop = (float)(i)/dataset.size();
        int progress = 51*prop;
        std::cout << "\r" << "[" << std::string(progress, '=') << std::string(50-progress, ' ') << "]  (" << prop*100.0 << "%)";
        std::cout.flush();

        #if NESTED_TYPE != 0
          } // end if (master thread)
        #endif
      }
      #if NESTED_TYPE != 0
         } // end omp parallel
         // Free local ged methods
         for  (unsigned int _g=0; _g<local_ged.size(); _g++)
           delete local_ged[_g];
      #endif

      //std::cout << "\r" << "[" << std::string(50, '=') << "]" << std::endl;
      std::cout << "\r  [" << GREEN << "done" << RESET << "]" << std::string(55, ' ') << std::endl;

      //update the weight of egdes with projection frequency
      //countFreqProj( sg, mappings, graphs);

      std::ofstream ofs_proj (out_proj, std::ofstream::out);
      outputJSON_proj(sg, mappings, graphs, id_graphs, dataset, &ofs_proj);
      ofs_proj.close();

      for (std::list<ECMapping*>::iterator it=mappings.begin(); it!=mappings.end(); it++)
        delete *it;

    } // End if(compute projections)


    /******************************************
     * Output the list of train, val and test indices
     ******************************************/

   if(!use_idFile){

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
   }


  } // end for each subset


  delete cf;
  delete ed_init;
  delete ipfp;
  delete ptr_dataset;
  //delete sg; // < deleted by the dataset !

  return 0;
}
