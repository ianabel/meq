

MFEM_DIR=/home/ian/projects/mfem/install

include $(MFEM_DIR)/share/mfem/config.mk

SUNDIALS_DIR=/home/ian/projects/sundials/install

MFEM_OPT = -I$(MFEM_DIR)/include $(MFEM_LIBS) -I$(SUNDIALS_DIR)/include -I/usr/include/suitesparse/

CXX = g++
CXXFLAGS_RELEASE = -std=c++17 -O3 -Wall -mtune=native
CXXFLAGS_DEBUG = -std=c++17 -g -Og -Wall

mfemPoisson: mfemPoisson.cpp PoissonSolver.cpp
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ mfemPoisson.cpp PoissonSolver.cpp $(MFEM_OPT)

clean:
	rm -f mfemPoisson

distclean: clean
	echo "Nothing to do"

.PHONY: clean distclean

