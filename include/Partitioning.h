#ifndef PARTITIONING_H
#define PARTITIONING_H

#include <metis.h>
#include <vector>
#include "Simulation.h"
#include <cstdlib>

class Partitioning {
public:
	idx_t* partitioning;
	
	Partitioning(std::vector<Population*>* pops, int population_size) {
		partitioning = (idx_t*)calloc(population_size,sizeof(idx_t));
		
	}
	virtual ~Partitioning() {
		free(partitioning);
	}
	
	/*virtual void clear_intermediate_data() {}	
	virtual void generate_connection_structures(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {}
	virtual void assign_partition_units(Model* model,int partitions, int process_id) {}	
	*/
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {}
	
};

#endif
