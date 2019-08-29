all: meq

MFEM_DIR=/home/ian/projects/ian-mfem
MFEM_LIB=-lmfem -lrt -llapack -lblas
MFEM_DEBUG_LIB=-lmfem-debug -lrt
SUNDIALS_LIB = -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol
CXX = g++
CXXFLAGS_RELEASE = -std=c++17 -O3 -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall

Baxis: Baxis.cpp
	$(CXX) $(CXXFLAGS_DEBUG) Baxis.cpp -o Baxis $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)

mfemSurfaceTracer: mfemSurfaceTracer.cpp
	$(CXX) $(CXXFLAGS_DEBUG) mfemSurfaceTracer.cpp -o mfemSurfaceTracer $(MFEM_DEBUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


