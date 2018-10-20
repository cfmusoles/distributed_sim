project := distSim

METIS_PATH=$(HOME)/metis_build
PARMETIS_PATH=$(HOME)/parmetis_ompi_build
ZOLTAN_PATH=$(HOME)/Zoltan_v3.83/build

CXX = mpicxx
SRC_PATH = ./src
SRC_EXT = cpp

FILES = $(shell find $(SRC_PATH) -name '*.$(SRC_EXT)' | sort -k 1nr |cut -f2-) 
OUT = $(project)

INCLUDES = -I./include -I$(PARMETIS_PATH)/include -I$(METIS_PATH)/include -I$(ZOLTAN_PATH)/include 
LIBS = -L$(PARMETIS_PATH)/lib -L$(METIS_PATH)/lib -L$(ZOLTAN_PATH)/lib
 
LDFLAGS = -lmetis -lparmetis -lzoltan
LSTATIC = 

build: $(FILES)
	$(CXX) -O3 -o $(OUT) $(INCLUDES) $(LIBS) $(FILES) $(LDFLAGS) $(LSTATIC) --std=c++11

.PHONY: debug
debug:
	$(CXX) -g -o $(OUT) $(INCLUDES) $(LIBS) $(FILES) $(LDFLAGS) $(LSTATIC) --std=c++11

.PHONY: clean
clean:
	rm $(OUT)

