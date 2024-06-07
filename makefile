CXX       = mpic++
CXXFLAGS ?= -std=c++20

#change this path
PACS_ROOT?=../pacs-examples/Examples

export DOXYFILE=./Doxyfile

CPPFLAGS ?= -fopenmp -O3 -Wall -Wno-conversion-null -Wno-deprecated-declarations -pedantic -I. -I$(PACS_ROOT)/src/Utilities/ -I$(PACS_ROOT)/include/  -I/mnt/wslg/distro/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3

LDFLAGS ?= -L$(PACS_ROOT)/src/Utilities/ -L$(PACS_ROOT)/lib/
LIBS  ?= -lpacs -lmuparser

# Get all files *.cpp
SRCS=$(wildcard *.cpp)
# Get the corresponding object file
OBJS = $(SRCS:.cpp=.o)
# Get all headers in the working directory
HEADERS=$(wildcard *.hpp)
#
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)


#========================== ORA LA DEFINIZIONE DEGLI OBIETTIVI
.PHONY: all clean distclean doc

.DEFAULT_GOAL = all

all: $(EXEC)

clean:
	-\rm -f $(EXEC) $(OBJS)

distclean: clean
	-\rm -f ./doc $(DEPEND)
	-\rm -f *.out *.bak *~

doc:
	doxygen $(DOXYFILE)

$(EXEC): $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

# Compile each source file into an object file
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@