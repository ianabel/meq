
all: meq

MFEM_DIR=/home/ian/projects/ian-mfem
MFEM_LIB=-lmfem -lrt -llapack -lblas
MFEM_DEBUG_LIB=-lmfem-debug -lrt
SUNDIALS_LIB = -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol
CXXFLAGS_RELEASE = -std=c++17 -fopenmp -Ofast -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall
GSINVERTER_DEP = GSInverter.cpp GSInverter.hpp HDGGSIntegrator.cpp HDGGSIntegrator.hpp CockburnEstimator.hpp
GSINVERTER_SRC = GSInverter.cpp HDGGSIntegrator.cpp



meq: meq.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) meq.cpp $(GSINVERTER_SRC) -o meq $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


meq-debug: meq.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_DEBUG) meq.cpp $(GSINVERTER_SRC) -o meq $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


mfemHDGPoisson: mfemHDGPoisson.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_DEBUG) mfemHDGPoisson.cpp $(GSINVERTER_SRC) -o mfemHDGPoisson $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemHDGPoissonNL: mfemHDGPoissonNL.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) mfemHDGPoissonNL.cpp $(GSINVERTER_SRC) -o mfemHDGPoissonNL $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemHDGPoissonNL-KINSOL: mfemHDGPoissonNL-KINSOL.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) mfemHDGPoissonNL-KINSOL.cpp $(GSINVERTER_SRC) -o mfemHDGPoissonNL-KINSOL $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfem-test: mfem-test.cpp
	g++ $(CXXFLAGS_DEBUG) mfem-test.cpp -o mfem-test $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


