
all: meq

MFEM_SOURCE_DIR=mfem

ifdef DEBUG
	MFEM_DIR = mfem-debug
else
	MFEM_DIR = mfem-release
endif

MFEM_CONFIG = $(MFEM_DIR)/config/config.mk

include $(MFEM_CONFIG)

SUNDIALS_DIR=$(CURDIR)/sundials-install

CXX = $(MFEM_CXX)
MEQ_FLAGS = -O3 -std=c++17 -Wall

GSINVERTER_DEP = GSSolver.cpp GSSolver.hpp HDGGSIntegrator.cpp HDGGSIntegrator.hpp CockburnEstimator.hpp FreeBoundary.hpp Solution.hpp
GSINVERTER_SRC = GSSolver.cpp HDGGSIntegrator.cpp

HEADER_DEP = meq.hpp config.hpp utility.hpp model.hpp

meq: meq.cpp FreeBoundary.cpp $(GSINVERTER_DEP) $(MFEM_LIB_FILE)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ meq.cpp FreeBoundary.cpp $(GSINVERTER_SRC) $(MFEM_LIBS)

mfemGS: mfemGS.cpp SolovievEquilibrium.hpp $(HEADER_DEP) $(GSINVERTER_DEP) $(MFEM_LIB_FILE)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ mfemGS.cpp $(GSINVERTER_SRC) $(MFEM_LIBS)

meshgenpp: meshgenpp.cpp meshgen.cpp $(HEADER_DEP) $(MFEM_LIB_FILE)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ meshgenpp.cpp meshgen.cpp -lgmsh $(MFEM_LIBS)

mfemProjector: mfemProjector.cpp $(HEADER_DEP) $(MFEM_LIB_FILE)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ mfemProjector.cpp  $(MFEM_LIBS) -lnetcdf_c++4 -lnetcdf

mfemCheck: mfemCheck.cpp $(HEADER_DEP) $(MFEM_LIB_FILE)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ mfemCheck.cpp $(MFEM_LIBS)

FluxSurfaces: FluxSurfaces.cpp $(HEADER_DEP) $(MFEM_LIB_FILE)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ FluxSurfaces.cpp $(MFEM_LIBS)


sundials/.git: 
	git submodule update --init sundials

SUNDIALS_BUILD_DIR = $(CURDIR)/sundials-build

SUNDIALS_CMAKE_OPT = -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX=$(SUNDIALS_DIR) -DEXAMPLES_INSTALL=off
$(SUNDIALS_BUILD_DIR)/Makefile: sundials/.git
	cmake $(SUNDIALS_CMAKE_OPT) -B $(SUNDIALS_BUILD_DIR) -S sundials

$(SUNDIALS_DIR)/include: $(SUNDIALS_BUILD_DIR)/Makefile
	rm -rf $(SUNDIALS_DIR); mkdir $(SUNDIALS_DIR); make $(MAKEFLAGS) -C $(SUNDIALS_BUILD_DIR) install;

toml11/toml.hpp:
	git submodule update --init toml11;

mfem/.git:
	git submodule update --init mfem;

ifdef DEBUG
MFEM_CMAKE = -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-g -Og" -DMFEM_USE_SUNDIALS=ON -DSUNDIALS_DIR=$(SUNDIALS_DIR)
else
MFEM_CMAKE = -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O3 -mtune=native -DNDEBUG" -DMFEM_USE_SUNDIALS=ON -DSUNDIALS_DIR=$(SUNDIALS_DIR)
endif

$(MFEM_CONFIG): mfem/.git $(SUNDIALS_DIR)/include
	cmake -B $(MFEM_DIR) -S$(MFEM_SOURCE_DIR) $(MFEM_CMAKE)

$(MFEM_LIB_FILE): $(MFEM_CONFIG)
	make $(MAKEFLAGS) -C $(MFEM_DIR)

vacuum-test: VacuumGFSoln.cpp FreeBoundary.cpp $(HEADER_DEP)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ VacuumGFSoln.cpp FreeBoundary.cpp $(MFEM_LIBS) -lnetcdf_c++4 -lnetcdf


vacuum-mesh-test: VacuumMeshSoln.cpp FreeBoundary.cpp $(HEADER_DEP)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ VacuumGFSoln.cpp FreeBoundary.cpp $(MFEM_LIBS)

clean:
	rm -f meshgenpp mfemGS mfemGS-debug meq meq-debug mfemCheck mfemProjector vacuum-test vacuum-mesh-test

distclean: clean
	rm -rf $(MFEM_DIR) $(SUNDIALS_BUILD_DIR) $(SUNDIALS_DIR)

.PHONY: clean distclean $(MFEM_LIB_FILE)

