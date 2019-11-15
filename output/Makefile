all: meq

MFEM_DIR=mfem
MFEM_RELEASE_DIR=../mfem-release
MFEM_DEBUG_DIR=../mfem-debug
MFEM_LIB= -L$(MFEM_RELEASE_DIR) -lmfem -lrt 
MFEM_DEBUG_LIB= -L$(MFEM_DEBUG_DIR) -lmfem -lrt

SUNDIALS_DIR=../sundials-install
SUNDIALS_LIB = -I$(SUNDIALS_DIR)/include -Wl,-rpath,$(SUNDIALS_DIR)/lib -L$(SUNDIALS_DIR)/lib -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol

MFEM_RELEASE_OPT = -I$(MFEM_RELEASE_DIR) $(MFEM_LIB) $(SUNDIALS_LIB)
MFEM_DEBUG_OPT   = -I$(MFEM_DEBUG_DIR) $(MFEM_LIB) $(SUNDIALS_LIB)

TOML_DIR=../toml11
TOML_OPT=-I$(TOML_DIR)

CXX = g++
CXXFLAGS_RELEASE = -std=c++17 -O3 -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall

Baxis: Baxis.cpp
	$(CXX) $(CXXFLAGS_DEBUG) Baxis.cpp -o Baxis $(MFEM_DEBUG_OPT) $(TOML_OPT)

mfemSurfaceTracer: mfemSurfaceTracer.cpp
	$(CXX) $(CXXFLAGS_DEBUG) mfemSurfaceTracer.cpp -o mfemSurfaceTracer $(MFEM_DEBUG_OPT)


