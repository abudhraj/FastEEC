CXX = g++
CXXFLAGS = -O3 -std=c++14 -pedantic -W -Wall -Wshadow -fPIC -pthread -ffast-math

# Find directory paths
FASTJETINC := $(shell fastjet-config --cxxflags)
FASTJETLIB := $(shell fastjet-config --libs)
INCLUDES += $(FASTJETINC) -Iinclude
LIBRARIES += $(FASTJETLIB) -ldl

SRCS := $(wildcard *.cc)
TARGETS := $(patsubst %.cc,%,$(SRCS))

all: $(TARGETS)

%: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(LIBRARIES)
	
clean:
	rm -f $(TARGETS) *.o
