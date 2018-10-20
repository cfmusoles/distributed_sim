#ifndef CUSTOM_PARTITIONING
#define CUSTOM_PARTITIONING

#include <metis.h>
#include <vector>
#include <algorithm>
#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"

class CustomPartitioning :  public Partitioning {
public:
	
	CustomPartitioning(std::vector<Population*>* pops, int population_size, idx_t* parts) : Partitioning(pops,population_size) {
        partitioning = parts;
    }
	virtual ~CustomPartitioning() {}
	
};

#endif