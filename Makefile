
all: meq

MFEM_DIR=/home/ian/projects/mfem-hdg-branch
MFEM_LIB=-lmfem -lrt
MFEM_LIB_DEBUG=-lmfem-debug -lrt
SUNDIALS_LIB = -lsundials_arkode -lsundials_cvode -lsundials_nvecserial -lsundials_kinsol
GS_INVERT_DEP = GSInverter.hpp GSInverter.cpp
GS_INVERT_SOURCE = GSInverter.cpp



meq: meq.cpp $(GS_INVERT_DEP)
	g++ -std=c++17 -O3 -Wall meq.cpp $(GS_INVERT_SOURCE) -o meq $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


meq-debug: meq.cpp
	g++ -std=c++17 -g -Og -Wall meq.cpp -o meq $(MFEM_BEDUG_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)


mfemHDGPoisson: mfemHDGPoisson.cpp $(GS_INVERT_DEP)
	g++ -std=c++17 -O3 -Wall mfemHDGPoisson.cpp $(GS_INVERT_DEP) -o mfemHDGPoisson $(MFEM_LIB) $(SUNDIALS_LIB) -I$(MFEM_DIR) -L$(MFEM_DIR)
