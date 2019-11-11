
all: meq

MFEM_DIR=./mfem
MFEM_RELEASE_DIR=./mfem-release
MFEM_DEBUG_DIR=./mfem-debug
MFEM_LIB= -L$(MFEM_RELEASE_DIR) -lmfem -lrt 
MFEM_DEBUG_LIB= -L$(MFEM_DEBUG_DIR) -lmfem -lrt

SUNDIALS_DIR=/home/ian/projects/sundials-install
SUNDIALS_LIB = -I$(SUNDIALS_DIR)/include -Wl,-rpath,$(SUNDIALS_DIR)/lib -L$(SUNDIALS_DIR)/lib -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol

MFEM_RELEASE_OPT = -I$(MFEM_RELEASE_DIR) $(MFEM_LIB) $(SUNDIALS_LIB)
MFEM_DEBUG_OPT   = -I$(MFEM_DEBUG_DIR) $(MFEM_LIB) $(SUNDIALS_LIB)

CXX = g++
CXXFLAGS_RELEASE = -std=c++17 -O3 -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall
GSINVERTER_DEP = GSSolver.cpp GSSolver.hpp HDGGSIntegrator.cpp HDGGSIntegrator.hpp CockburnEstimator.hpp FreeBoundary.hpp Solution.hpp
GSINVERTER_SRC = GSSolver.cpp HDGGSIntegrator.cpp

MFEM_RELEASE_DEP = $(MFEM_RELEASE_DIR)/libmfem.a
MFEM_DEBUG_DEP = $(MFEM_DEBUG_DIR)/libmfem.a

MAKEFLAGS ?= -j4

meq: meq.cpp FreeBoundary.cpp $(GSINVERTER_DEP) $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ meq.cpp FreeBoundary.cpp $(GSINVERTER_SRC) $(MFEM_RELEASE_OPT)

meq-debug: meq.cpp $(GSINVERTER_DEP) $(MFEM_DEBUG_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ meq.cpp FreeBoundary.cpp $(GSINVERTER_SRC) $(MFEM_DEBUG_OPT)

mfemGS: mfemGS.cpp StdFnCoeffs.hpp SolovievEquilibrium.hpp $(GSINVERTER_DEP) $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ mfemGS.cpp $(GSINVERTER_SRC) $(MFEM_RELEASE_OPT)

mfemGS-debug: mfemGS.cpp StdFnCoeffs.hpp SolovievEquilibrium.hpp $(GSINVERTER_DEP) $(MFEM_DEBUG_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ mfemGS.cpp $(GSINVERTER_SRC) $(MFEM_DEBUG_OPT)

meshgenpp: meshgenpp.cpp meshgen.cpp MEQConf.hpp meq.hpp $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ meshgenpp.cpp meshgen.cpp -lgmsh $(MFEM_RELEASE_DEP)

mfemProjector: mfemProjector.cpp MEQConf.hpp meq.hpp $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ mfemProjector.cpp  $(MFEM_RELEASE_OPT)

mfemCheck: mfemCheck.cpp MEQConf.hpp meq.hpp $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ mfemCheck.cpp $(MFEM_RELEASE_OPT)

toml11/toml.hpp:
	git submodule init toml11\
	git submodule update toml11

mfem/.git:
	git submodule init mfem\
	git submodule update mfem

MFEM_RELEASE_CMAKE = -DCMAKE_BUILD_TYPE=Release -DMFEM_USE_SUNDIALS=ON -DSUNDIALS_DIR=$(SUNDIALS_DIR)
MFEM_DEBUG_CMAKE = -DCMAKE_BUILD_TYPE=Release -DMFEM_USE_SUNDIALS=ON -DSUNDIALS_DIR=$(SUNDIALS_DIR)

$(MFEM_RELEASE_DIR)/Makefile: |mfem/.git
	cmake -B $(MFEM_RELEASE_DIR) -S$(MFEM_DIR) $(MFEM_RELEASE_CMAKE)

mfem-debug/Makefile: |mfem/.git
	cmake -B $(MFEM_DEBUG_DIR) -S$(MFEM_DIR) $(MFEM_DEBUG_CMAKE)

$(MFEM_RELEASE_DIR)/libmfem.a: $(MFEM_RELEASE_DIR)/Makefile
	make -C $(MFEM_RELEASE_DIR)

$(MFEM_DEBUG_DIR)/libmfem.a: $(MFEM_DEBUG_DIR)/Makefile
	make -C $(MFEM_DEBUG_DIR)

clean:
	rm -f meshgenpp mfemGS mfemGS-debug meq meq-debug mfemCheck mfemProjector

distclean: clean
	rm -f $(MFEM_RELEASE_DIR) $(MFEM_DEBUG_DIR)

.PHONY: clean distclean $(MFEM_RELEASE_DIR)/libmfem.a $(MFEM_DEBUG_DIR)/libmfem.a

