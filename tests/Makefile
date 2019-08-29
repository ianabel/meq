# 
# Copyright Ian G Abel 2018
#

# If boost is installed site-wide, then don't set anything
# If you have installed a local boost install, set BOOSTDIR to the location
# of the install, or you can set BOOSTINCDIR and BOOSTLIBDIR manually

ifdef BOOSTDIR
BOOSTINCDIR = $(BOOSTDIR)/include
BOOSTLIBDIR = $(BOOSTDIR)/lib
endif

ifdef BOOSTINCDIR
CXXFLAGS += -I$(BOOSTINCDIR) -L$(BOOSTLIBDIR) -Wl,-rpath-link=$(BOOSTLIBDIR)  
endif

CXXFLAGS += -Wall -pedantic -std=c++17
CXXDEBUGFLAGS +=  -O0 -g
CXXTESTFLAGS +=   -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework
CXXRELEASEFLAGS += -O3


all: test

test: AnalyticGSTests FreeBoundaryTests ConfigTests HDGTests
	LD_LIBRARY_PATH=$(BOOSTLIBDIR) ./AnalyticGSTests
	LD_LIBRARY_PATH=$(BOOSTLIBDIR) ./FreeBoundaryTests
	LD_LIBRARY_PATH=$(BOOSTLIBDIR) ./ConfigTests
	LD_LIBRARY_PATH=$(BOOSTLIBDIR) ./HDGTests


AnalyticGSTests: AnalyticGSTests.cpp 
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXTESTFLAGS) -o $@ AnalyticGSTests.cpp

FreeBoundaryTests: FreeBoundaryTests.cpp 
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXTESTFLAGS) -o $@ FreeBoundaryTests.cpp

ConfigTests: ConfigTests.cpp 
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXTESTFLAGS) -o $@ ConfigTests.cpp
	
HDGTests: HDGTests.cpp 
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAGS) $(CXXTESTFLAGS) -o $@ HDGTests.cpp


clean:
	rm -f AnalyticGSTests FreeBoundaryTests ConfigTests HDGTests

.PHONY: clean test
