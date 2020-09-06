LSAPE_DIR=~/dev/lsape/include/
GRAPH_DIR=~/dev/graph-lib/include/
GRAPH_OBJ=~/dev/graph-lib/obj/
EIGEN_DIR=/usr/include/eigen3/
SINKHORN_DIR=~/dev/lsape_sinkhorn/
EDMONDS_DIR=~/dev/edmonds/
EDMONDS_LIB=edmonds.a

CXXFLAGS = -I$(LSAPE_DIR) -I$(GRAPH_DIR) -I$(EIGEN_DIR) -I$(EDMONDS_DIR) -I$(SINKHORN_DIR) -Wall  -std=c++11 -g -fopenmp
OBJ = $(GRAPH_OBJ)utils.o $(GRAPH_OBJ)ConstantGraphEditDistance.o $(GRAPH_OBJ)SymbolicGraph.o edmond.a

SRC=src
LIB=lib
INCL=include
BINDIR=bin
TESTDIR=tests


# Number of threads by default (should be <= number of processors)
ifeq ($(THREADS), "")
  THREADS=40
endif

$(TESTDIR)/test_supergraph: $(SRC)/test-supergraph.cpp $(INCL)/supergraph.h $(INCL)/utils.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -ltinyxml -I$(INCL)

$(BINDIR)/sg_json: $(SRC)/mk-supergraph-json.cpp $(INCL)/supergraph.h $(INCL)/SupergraphCostFunction.h $(INCL)/utils.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -ltinyxml -I$(INCL)


supergraph-cv-bipartite: CXXFLAGS += -DBIPARTITE_MATRIX
supergraph-cv-bipartite: $(BINDIR)/supergraph-cv-omp


$(BINDIR)/supergraph: $(SRC)/compute-supergraph.cpp $(INCL)/supergraph.h $(INCL)/dataset_json.h $(INCL)/utils.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -ltinyxml -I. -I../rapidjson/include

$(BINDIR)/supergraph-cv: $(SRC)/compute-supergraph-cv.cpp $(INCL)/supergraph.h $(INCL)/dataset_json.h $(INCL)/utils.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -ltinyxml -I. -I../rapidjson/include

$(BINDIR)/supergraph-cv-omp: $(SRC)/compute-supergraph-cv.cpp $(INCL)/supergraph.h $(INCL)/dataset_json.h $(INCL)/utils.h
	$(CXX) -o $@ $< $(CXXFLAGS) -DNESTED_TYPE=0 -DNUM_THREADS=$(THREADS) $(OBJ) -ltinyxml -I. -I../rapidjson/include



$(BINDIR)/project: $(SRC)/project.cpp $(INCL)/supergraph.h  $(INCL)/dataset_json.h $(INCL)/utils.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -I../rapidjson/include -I$(INCL) -ltinyxml

$(BINDIR)/stars: $(SRC)/compute-stars.cpp $(INCL)/utils.h $(INCL)/stars.h $(INCL)/dataset_json.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -ltinyxml -I. -I../rapidjson/include

$(BINDIR)/ds2json: utils/dataset2json.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -I$(INCL) -ltinyxml -I.

$(BINDIR)/ds_stats: utils/dataset_stats.cpp $(INCL)/utils.h $(INCL)/dataset_json.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -ltinyxml -I. -I../rapidjson/include

$(BINDIR)/cv-indices: utils/cv-indices.cpp $(INCL)/utils.h $(INCL)/dataset_json.h
	$(CXX) -o $@ $< $(CXXFLAGS) $(OBJ) -ltinyxml -I. -I../rapidjson/include


$(EDMONDS): $(EDMONDS_DIR)EdmondsMatching.o $(EDMONDS_DIR)Blossom.o $(EDMONDS_DIR)Vertex.o  $(EDMONDS_DIR)Edge.o
	ar rcs $@ $^

$(EDMONDS)%.o: $(EDMONDS_DIR)%.cpp $(EDMONDS_DIR)%.h
	$(CXX) -c $< -o $@ -I$(EDMONDS_DIR)


