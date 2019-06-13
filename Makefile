
all: meq

MFEM_DIR=/home/ian/projects/mfem-hdg-branch
MFEM_LIB=-lmfem -lrt
MFEM_LIB_DEBUG=-lmfem-debug -lrt
SUNDIALS_LIB = -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol
CXXFLAGS_RELEASE = -std=c++17 -Ofast -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall
GSINVERTER_DEP = GSInverter.cpp GSInverter.hpp HDGGSIntegrator.cpp HDGGSIntegrator.hpp
GSINVERTER_SRC = GSInverter.cpp HDGGSIntegrator.cpp



meq: meq.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) meq.cpp $(GSINVERTER_SRC) -o meq $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


meq-debug: meq.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_DEBUG) meq.cpp $(GSINVERTER_SRC) -o meq $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


mfemHDGPoisson: mfemHDGPoisson.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) mfemHDGPoisson.cpp $(GSINVERTER_SRC) -o mfemHDGPoisson $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemHDGPoissonNL: mfemHDGPoissonNL.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) mfemHDGPoissonNL.cpp $(GSINVERTER_SRC) -o mfemHDGPoissonNL $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)
