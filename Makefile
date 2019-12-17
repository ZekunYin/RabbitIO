DIR_INC := .
DIR_SRC := .
DIR_OBJ := .

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := fastIO

BIN_TARGET := ${TARGET}

CXX := g++
CXXFLAGS := -std=c++11 -g -fopenmp -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LIBS := -lz -lpthread -fopenmp 
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)


${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)
	cp $(BIN_TARGET) ../

%.o:%.cpp 
	$(CXX) $(CXXFLAGS) -O3 -c $< -o $@

.PHONY:clean
clean:
	rm *.o $(TARGET) ../$(TARGET)

#make_obj_dir:
#	@if test ! -d $(DIR_OBJ) ; \
#	then \
#		mkdir $(DIR_OBJ) ; \
#	fi

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."
