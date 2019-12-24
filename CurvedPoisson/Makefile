

MFEM_DIR=../mfem-perf
MFEM_CONFIG = $(MFEM_DIR)/config/config.mk

include $(MFEM_CONFIG)


CXX = $(MFEM_CXX)
MEQ_FLAGS = -std=c++17 -Wall


mfemPoisson: mfemPoisson.cpp PoissonSolver.cpp $(MFEM_LIB_FILE)
	$(CXX) $(MEQ_FLAGS) $(MFEM_FLAGS) -o $@ mfemPoisson.cpp PoissonSolver.cpp $(MFEM_LIBS)

clean:
	rm -f mfemPoisson

distclean: clean
	echo "Nothing to do"

.PHONY: clean distclean 

