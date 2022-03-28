
all: meq

include Makefile.local

MFEM_CONFIG = $(MFEM_INSTALL_DIR)/share/mfem/config.mk

include $(MFEM_CONFIG)

CXX = $(MFEM_CXX)
MEQ_FLAGS = -O3 -std=c++17 -Wall -I$(TOML11_DIR)

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
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ FluxSurfaces.cpp $(MFEM_LIBS) -lnetcdf_c++4 -lnetcdf

vacuum-test: VacuumGFSoln.cpp FreeBoundary.cpp $(HEADER_DEP)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ VacuumGFSoln.cpp FreeBoundary.cpp $(MFEM_LIBS) -lnetcdf_c++4 -lnetcdf

vacuum-mesh-test: VacuumMeshSoln.cpp FreeBoundary.cpp $(HEADER_DEP)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ VacuumGFSoln.cpp FreeBoundary.cpp $(MFEM_LIBS)

clean:
	rm -f meshgenpp mfemGS mfemGS-debug meq meq-debug mfemCheck mfemProjector vacuum-test vacuum-mesh-test

distclean: clean
	rm -rf $(MFEM_DIR) $(SUNDIALS_BUILD_DIR) $(SUNDIALS_DIR)

.PHONY: clean distclean $(MFEM_LIB_FILE)

