

MFEM_RELEASE_DIR=../mfem-release
MFEM_DEBUG_DIR=../mfem-debug
MFEM_LIB= -L$(MFEM_RELEASE_DIR) -lmfem -lrt 
MFEM_DEBUG_LIB= -L$(MFEM_DEBUG_DIR) -lmfem -lrt

SUNDIALS_DIR=../sundials-install
SUNDIALS_LIB = -I$(SUNDIALS_DIR)/include -Wl,-rpath,$(SUNDIALS_DIR)/lib -L$(SUNDIALS_DIR)/lib -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol

MFEM_RELEASE_OPT = -I$(MFEM_RELEASE_DIR) $(MFEM_LIB) $(SUNDIALS_LIB)
MFEM_DEBUG_OPT   = -I$(MFEM_DEBUG_DIR) $(MFEM_LIB) $(SUNDIALS_LIB)

CXX = g++
CXXFLAGS_RELEASE = -std=c++17 -O3 -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall

MFEM_RELEASE_DEP = $(MFEM_RELEASE_DIR)/libmfem.a
MFEM_DEBUG_DEP = $(MFEM_DEBUG_DIR)/libmfem.a

mfemPoisson: mfemPoisson.cpp PoissonSolver.cpp
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ mfemPoisson.cpp PoissonSolver.cpp $(MFEM_RELEASE_OPT)

mfemPoisson-prof: mfemPoisson.cpp PoissonSolver.cpp
	$(CXX) $(CXXFLAGS_RELEASE) -pg -o $@ mfemPoisson.cpp PoissonSolver.cpp $(MFEM_RELEASE_OPT)


mfemPoisson-debug: mfemPoisson.cpp PoissonSolver.cpp
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ mfemPoisson.cpp PoissonSolver.cpp $(MFEM_DEBUG_OPT)

../mfem-release/libmfem.a:
	make $(MAKEFLAGS) -C ../ mfem-release/libmfem.a

../mfem-debug/libmfem.a:
	make $(MAKEFLAGS) -C ../ mfem-debug/libmfem.a

clean:
	rm -f mfemPoisson

distclean: clean
	echo "Nothing to do"

.PHONY: clean distclean ../mfem-release/libmfem.a ../mfem-debug/libmfem.a

