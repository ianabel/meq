
all: meq

MFEM_DIR=mfem
MFEM_RELEASE_DIR=./mfem-release
MFEM_DEBUG_DIR=./mfem-debug
MFEM_LIB= -L$(MFEM_RELEASE_DIR) -lmfem -lrt 
MFEM_DEBUG_LIB= -L$(MFEM_DEBUG_DIR) -lmfem -lrt

SUNDIALS_DIR=$(CURDIR)/sundials-install
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

HEADER_DEP = meq.hpp config.hpp utility.hpp model.hpp

MAKEFLAGS ?= -j4

meq: meq.cpp FreeBoundary.cpp $(GSINVERTER_DEP) $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ meq.cpp FreeBoundary.cpp $(GSINVERTER_SRC) $(MFEM_RELEASE_OPT)

meq-debug: meq.cpp $(GSINVERTER_DEP) $(MFEM_DEBUG_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ meq.cpp FreeBoundary.cpp $(GSINVERTER_SRC) $(MFEM_DEBUG_OPT)

mfemGS: mfemGS.cpp SolovievEquilibrium.hpp $(HEADER_DEP) $(GSINVERTER_DEP) $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ mfemGS.cpp $(GSINVERTER_SRC) $(MFEM_RELEASE_OPT)

mfemGS-debug: mfemGS.cpp SolovievEquilibrium.hpp $(HEADER_DEP) $(GSINVERTER_DEP) $(MFEM_DEBUG_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ mfemGS.cpp $(GSINVERTER_SRC) $(MFEM_DEBUG_OPT)

meshgenpp: meshgenpp.cpp meshgen.cpp $(HEADER_DEP) $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ meshgenpp.cpp meshgen.cpp -lgmsh $(MFEM_RELEASE_OPT)

mfemProjector: mfemProjector.cpp $(HEADER_DEP) $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ mfemProjector.cpp  $(MFEM_RELEASE_OPT)

mfemCheck: mfemCheck.cpp $(HEADER_DEP) $(MFEM_RELEASE_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ mfemCheck.cpp $(MFEM_RELEASE_OPT)

sundials/.git: 
	git submodule update --init sundials

SUNDIALS_BUILD_DIR = $(CURDIR)/sundials-build

SUNDIALS_CMAKE_OPT = -DCMAKE_INSTALL_PREFIX=$(SUNDIALS_DIR) -DEXAMPLES_INSTALL=off
$(SUNDIALS_BUILD_DIR)/Makefile: sundials/.git
	cmake $(SUNDIALS_CMAKE_OPT) -B $(SUNDIALS_BUILD_DIR) -S sundials

$(SUNDIALS_DIR)/include: $(SUNDIALS_BUILD_DIR)/Makefile
	rm -rf $(SUNDIALS_DIR); mkdir $(SUNDIALS_DIR); make $(MAKEFLAGS) -C $(SUNDIALS_BUILD_DIR) install;

toml11/toml.hpp:
	git submodule update --init toml11;

mfem/.git:
	git submodule update --init mfem;

MFEM_RELEASE_CMAKE = -DCMAKE_BUILD_TYPE=Release -DMFEM_USE_SUNDIALS=ON -DSUNDIALS_DIR=$(SUNDIALS_DIR)
MFEM_DEBUG_CMAKE = -DCMAKE_BUILD_TYPE=Release -DMFEM_USE_SUNDIALS=ON -DSUNDIALS_DIR=$(SUNDIALS_DIR)

mfem-release/Makefile: mfem/.git $(SUNDIALS_DIR)/include
	cmake -B $(MFEM_RELEASE_DIR) -S$(MFEM_DIR) $(MFEM_RELEASE_CMAKE)

mfem-debug/Makefile: mfem/.git $(SUNDIALS_DIR)/include
	cmake -B $(MFEM_DEBUG_DIR) -S$(MFEM_DIR) $(MFEM_DEBUG_CMAKE)

mfem-release/libmfem.a: mfem-release/Makefile $(SUNDIALS_DIR)/include
	make $(MAKEFLAGS) -C mfem-release

mfem-debug/libmfem.a: mfem-debug/Makefile $(SUNDIALS_DIR)/include
	make $(MAKEFLAGS) -C mfem-debug

vacuum-test: VacuumGFSoln.cpp FreeBoundary.cpp $(HEADER_DEP) $(MFEM_DEBUG_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ VacuumGFSoln.cpp FreeBoundary.cpp $(MFEM_DEBUG_OPT)

vacuum-mesh-test: VacuumMeshSoln.cpp FreeBoundary.cpp $(HEADER_DEP) $(MFEM_DEBUG_DEP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ VacuumGFSoln.cpp FreeBoundary.cpp $(MFEM_DEBUG_OPT)


clean:
	rm -f meshgenpp mfemGS mfemGS-debug meq meq-debug mfemCheck mfemProjector

distclean: clean
	rm -rf $(MFEM_RELEASE_DIR) $(MFEM_DEBUG_DIR) $(SUNDIALS_BUILD_DIR) $(SUNDIALS_DIR)

.PHONY: clean distclean mfem-release/libmfem.a mfem-debug/libmfem.a

