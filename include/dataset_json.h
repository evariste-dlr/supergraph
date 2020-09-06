#ifndef __DATASET_JSON_H__
#define __DATASET_JSON_H__

#include <iostream>
#include <cstdio>

#include "Dataset.h"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include "include/utils.h"

/**
 * @brief  Fill in the given dataset with the data in file
 */
void dataset_json(char* file, ChemicalDataset<double> & dataset)
{
  FILE * fp = fopen(file, "rb");
  char readBuffer[65536];
  rapidjson::FileReadStream ifs(fp, readBuffer, sizeof(readBuffer));
  rapidjson::Document d;
  d.ParseStream(ifs);
  if (d.HasParseError()){
    std::cout << RED << "[E]  Parse error in the JSON file" << RESET << std::endl;
    fclose(fp);
    exit(1);
  }

  assert(d.IsArray());

  for (rapidjson::Value::ConstValueIterator itr = d.Begin(); itr != d.End(); ++itr){
    const rapidjson::Value& j_nodes = (*itr)["nodes"];
    const rapidjson::Value& j_adj = (*itr)["adj"];
    double value = (*itr)["value"].GetDouble();
    assert(j_nodes.IsArray());
    assert(j_adj.IsArray());

    Graph<int,int> * g = new Graph<int,int>();
    // Nodes
    int j=0;
    for (rapidjson::SizeType rj=0; rj<j_nodes.Size(); rj++, j++){
      const rapidjson::Value& val = j_nodes[rj];
      int n = val.GetInt();
      g->Add(new GNode<int,int>(j, n));
    }

    // Edges
    j=0;
    for (rapidjson::SizeType rj=0; rj<j_adj.Size(); rj++, j++){
      int k=j;
      for (rapidjson::SizeType rk=rj; rk<j_nodes.Size(); rk++, k++){
        const rapidjson::Value& val = j_adj[rj][rk];
        int a = val.GetInt();
        if (a != 0)
          g->Link(j,k,a);
      }
    }

    dataset.add(g, value);
  }

  fclose(fp);
}




/**
 * @brief  Load the indices of the json file `file_name` in the given vectors
 */
bool load_subsets( const std::string& file_name,
                   std::vector<uint>& trainIdx,
                   std::vector<uint>& testIdx,
                   std::vector<uint>& valIdx )
{
  FILE * fp = fopen(file_name.c_str(), "rb");
  char readBuffer[65536];
  rapidjson::FileReadStream ifs(fp, readBuffer, sizeof(readBuffer));
  rapidjson::Document d;
  d.ParseStream(ifs);
  if (d.HasParseError()){
    std::cerr << RED << "[E]  Parse error in the JSON file" << RESET << std::endl;
    fclose(fp);
    return false;
  }

  try{

    assert(d.IsObject());
    assert(d.HasMember("train"));
    assert(d.HasMember("val"));
    assert(d.HasMember("test"));

    const rapidjson::Value& train = d["train"];
    assert(train.IsArray());
    trainIdx.resize(train.Size());
    for (rapidjson::SizeType i = 0; i < train.Size(); i++)
      trainIdx[i] = train[i].GetInt();

    const rapidjson::Value& val = d["val"];
    assert(val.IsArray());
    valIdx.resize(val.Size());
    for (rapidjson::SizeType i = 0; i < val.Size(); i++)
      valIdx[i] = val[i].GetInt();

    const rapidjson::Value& test = d["test"];
    assert(test.IsArray());
    testIdx.resize(test.Size());
    for (rapidjson::SizeType i = 0; i < test.Size(); i++)
      testIdx[i] = test[i].GetInt();

  }
  catch (std::exception& e){
    std::cerr << RED << "[E]  Unable to read JSON file" << RESET << std::endl;
  }
  return true;
}


#endif //__DATASET_JSON_H__





