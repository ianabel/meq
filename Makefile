
all: meq

MFEM_DIR=/home/ian/projects/ian-mfem
MFEM_LIB=-lmfem -lrt -llapack -lblas
MFEM_DEBUG_LIB=-lmfem-debug -lrt
SUNDIALS_LIB = -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol
CXXFLAGS_RELEASE = -std=c++17 -O3 -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall
GSINVERTER_DEP = GSInverter.cpp GSInverter.hpp HDGGSIntegrator.cpp HDGGSIntegrator.hpp CockburnEstimator.hpp FreeBoundary.hpp
GSINVERTER_SRC = GSInverter.cpp HDGGSIntegrator.cpp



meq: meq.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) meq.cpp $(GSINVERTER_SRC) -o meq $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


meq-debug: meq.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_DEBUG) meq.cpp $(GSINVERTER_SRC) -o meq-debug $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


mfemHDGPoisson: mfemHDGPoisson.cpp $(GSINVERTER_DEP) FreeBoundary.hpp 
	g++ $(CXXFLAGS_DEBUG) mfemHDGPoisson.cpp $(GSINVERTER_SRC) -o mfemHDGPoisson $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemHDGPoissonNL: mfemHDGPoissonNL.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) mfemHDGPoissonNL.cpp $(GSINVERTER_SRC) -o mfemHDGPoissonNL $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemHDGPoissonNL-KINSOL: mfemHDGPoissonNL-KINSOL.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) mfemHDGPoissonNL-KINSOL.cpp $(GSINVERTER_SRC) -o mfemHDGPoissonNL-KINSOL $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemHDGPoissonNL-KINSOL-DEBUG: mfemHDGPoissonNL-KINSOL.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_DEBUG) mfemHDGPoissonNL-KINSOL.cpp $(GSINVERTER_SRC) -o mfemHDGPoissonNL-KINSOL-DEBUG $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemGS: mfemGS.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_RELEASE) mfemGS.cpp $(GSINVERTER_SRC) -o mfemGS $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemGS-DEBUG: mfemGS.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_DEBUG) mfemGS.cpp $(GSINVERTER_SRC) -o mfemGS-DEBUG $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)



mfem-test: mfem-test.cpp $(GSINVERTER_DEP)
	g++ $(CXXFLAGS_DEBUG) mfem-test.cpp HDGGSIntegrator.cpp -o mfem-test $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemSurfaceTracer: mfemSurfaceTracer.cpp
	g++ $(CXXFLAGS_DEBUG) mfemSurfaceTracer.cpp -o mfemSurfaceTracer $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

Baxis: Baxis.cpp
	g++ $(CXXFLAGS_DEBUG) Baxis.cpp -o Baxis $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

