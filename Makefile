project := distSim

PREP = scorep --user #instrumentation
CXX = mpicxx
ARCHER_CC = CC
SRC_PATH = ./src
SRC_EXT = cpp

FILES = $(shell find $(SRC_PATH) -name '*.$(SRC_EXT)' | sort -k 1nr |cut -f2-) 
OUT = $(project)

INCLUDES = -I./include -I../PhD/PRAW/include -I$(HOME)/parmetis_ompi_build/include -I$(HOME)/metis_build/include -I$(HOME)/scorep-3.0/include -I$(HOME)/Zoltan_v3.83/build/include 
LIBS = -L$(HOME)/parmetis_ompi_build/lib -L$(HOME)/metis_build/lib -L$(HOME)/Zoltan_v3.83/build/lib
 
INCLUDES_INTEL = -I./include -I../PhD/PRAW/include -I$(HOME)/parmetis_intel_build/include -I$(HOME)/metis_build/include -I$(HOME)/scorep-3.0/include -I$(HOME)/Zoltan_v3.83_intel/build/include 
INTEL_LIBS = -L$(HOME)/metis_build/lib -L$(HOME)/parmetis_intel_build/lib  -L$(HOME)/Zoltan_v3.83_intel/build/lib 

INCLUDES_ARCHER = -I./include -I../PhD/PRAW/include -I$(HOME)/parmetis-4.0.3/build/include -I$(HOME)/metis-5.1.0/build/include -I$(HOME)/Zoltan_v3.83/build/include 
ARCHER_LIBS = -L$(HOME)/metis-5.1.0/build/lib -L$(HOME)/parmetis-4.0.3/build/lib -L$(HOME)/Zoltan_v3.83/build/lib 

LDFLAGS = -lmetis -lparmetis -lzoltan
LSTATIC = 

build: $(FILES)
	$(CXX) -O3 -o $(OUT) $(INCLUDES) $(LIBS) $(FILES) $(LDFLAGS) $(LSTATIC) --std=c++11

.PHONY: intelmpi
intelmpi:
	bash intelmpisourcewrap.sh $(CXX) -O3 -o $(OUT) $(INCLUDES_INTEL) $(INTEL_LIBS) $(FILES) $(LDFLAGS) $(LSTATIC) --std=c++11

.PHONY: debug
debug:
	$(CXX) -g -o $(OUT) $(INCLUDES) $(LIBS) $(FILES) $(LDFLAGS) $(LSTATIC) --std=c++11

.PHONY: profile
profile:
	$(PREP) $(CXX) -g -o $(OUT) $(INCLUDES) $(LIBS) $(FILES) $(LDFLAGS) $(LSTATIC)

.PHONY: archer
archer:
	$(ARCHER_CC) -O3 -o $(OUT) $(INCLUDES_ARCHER) $(ARCHER_LIBS) $(FILES) $(LDFLAGS) $(LSTATIC)


.PHONY: clean
clean:
	rm $(OUT)

