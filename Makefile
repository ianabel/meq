
all: meq

MFEM_DIR=/home/ian/projects/mfem-hdg-ian
MFEM_LIB=-lmfem -lrt 
MFEM_DEBUG_LIB=-lmfem-debug -lrt
SUNDIALS_LIB = -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol
CXX = g++
CXXFLAGS_RELEASE = -std=c++17 -O3 -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall
GSINVERTER_DEP = GSInverter.cpp GSInverter.hpp HDGGSIntegrator.cpp HDGGSIntegrator.hpp CockburnEstimator.hpp FreeBoundary.hpp
GSINVERTER_SRC = GSInverter.cpp HDGGSIntegrator.cpp

meq: meq.cpp $(GSINVERTER_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) meq.cpp $(GSINVERTER_SRC) -o meq $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

meq-debug: meq.cpp $(GSINVERTER_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) meq.cpp $(GSINVERTER_SRC) -o meq-debug $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemGS: mfemGS.cpp StdFnCoeffs.hpp SolovievEquilibrium.hpp $(GSINVERTER_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) mfemGS.cpp $(GSINVERTER_SRC) -o mfemGS $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemGS-DEBUG: mfemGS.cpp StdFnCoeffs.hpp SolovievEquilibrium.hpp  $(GSINVERTER_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) mfemGS.cpp $(GSINVERTER_SRC) -o mfemGS-DEBUG $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


meshgenpp: meshgenpp.cpp meshgen.cpp MEQConf.hpp meq.hpp
	$(CXX) $(CXXFLAGS_DEBUG) meshgenpp.cpp meshgen.cpp -o meshgenpp -lgmsh $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

clean:
	rm -f meshgenpp mfemGS mfemGS-DEBUG meq meq-debug

.PHONY: clean

