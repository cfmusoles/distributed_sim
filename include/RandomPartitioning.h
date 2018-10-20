#ifndef RANDOM_PARTITIONING
#define RANDOM_PARTITIONING

#include <metis.h>
#include <vector>
#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"

class RandomPartitioning :  public Partitioning {
public:
	
	RandomPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {}
	virtual ~RandomPartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		idx_t nParts = partitions; // number of processes (partitions)
		
		if(nParts == 1) {
			PRINTF("Partitioning not required\n");
			for(int ii=0; ii < model->population_size; ii++) 
				partitioning[ii] = 0;
		} else {
			// randomise partition assignment (test)
			PRINTF("Random partitioning...\n");
			for(int ii=0; ii < model->population_size; ii++) {
				int r = floor(rand01() * (nParts));
				if(r >= nParts) r = nParts - 1;
				partitioning[ii] = r;
			}
		}
	}
};


/*class RandomPartitioning :  public Partitioning {
public:
	
	RandomPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {}
	virtual ~RandomPartitioning() {}
	
	virtual void assign_partition_units(Model* model, int partitions, int process_id) {
		idx_t nParts = partitions; // number of processes (partitions)
		
		if(nParts == 1) {
			PRINTF("Partitioning not required\n");
			for(int ii=0; ii < model->population_size; ii++) 
				partitioning[ii] = 0;
		} else {
			// randomise partition assignment (test)
			PRINTF("Random partitioning...\n");
			for(int ii=0; ii < model->population_size; ii++) {
				int r = floor(rand01() * (nParts));
				if(r >= nParts) r = nParts - 1;
				partitioning[ii] = r;
			}
		}
	}
};*/

#endif
