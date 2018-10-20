#ifndef ROUNDROBIN_PARTITIONING
#define ROUNDROBIN_PARTITIONING

#include <metis.h>
#include <vector>
#include <cmath>
#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"

class RoundRobinPartitioning :  public Partitioning {
public:
	
	RoundRobinPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {
        populations = pops;
    }
	virtual ~RoundRobinPartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		idx_t nParts = partitions; // number of processes (partitions)
		
		if(nParts == 1) {
			PRINTF("Partitioning not required\n");
			for(int ii=0; ii < model->population_size; ii++) 
				partitioning[ii] = 0;
		} else {
			// round robin assignment, based on population belonging
			PRINTF("Round robin partitioning...\n");
            
			// generate datastructure holding info of #incoming synapses for each neuron
			int* incoming_connections = (int*)calloc(model->population_size,sizeof(int));
			for(int ii=0; ii < model->population_size; ii++) {
				for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
					incoming_connections[abs(model->interconnections[ii][jj])] += 1;
				}
			}

			// calculate expected weight per process (number of syns + neurons)
			float expected_weight = ceil((model->total_connections + model->population_size) / partitions);
			
			// ROUND ROBIN ALLOCATION
			int current_partition = 0;
			int* actual_weights = (int*)calloc(partitions,sizeof(int));
			for(int ii=0; ii < populations->size(); ii++) {
                for(int nn=0; nn < populations->at(ii)->neuron_ids.size(); nn++) {
					int neuronid = populations->at(ii)->neuron_ids[nn];
                    partitioning[neuronid] = current_partition;
					actual_weights[current_partition] += incoming_connections[neuronid] + 1;
					do {
						current_partition++;
						if(current_partition >= partitions) current_partition = 0;
					} while(actual_weights[current_partition] >= expected_weight);
					
                }
			}
			free(actual_weights);

			

			/*
			// This method fills up a node before going to the next
			// for each partition
			//	pick next neuron from population and add it to the partition
			//	substract its weight (number of synapses) from the expected weight 
			//	do until actual weight <= expected weight
			int current_partition = 0;
            float actual_weight = 0;
			for(int ii=0; ii < populations->size(); ii++) {
                for(int nn=0; nn < populations->at(ii)->neuron_ids.size(); nn++) {
					int neuronid = populations->at(ii)->neuron_ids[nn];
                    if(current_partition < partitions-1 && actual_weight + incoming_connections[neuronid] * 0.5f > expected_weight) {
                        current_partition++;
						actual_weight = 0;
                    }
					partitioning[neuronid] = current_partition;
					actual_weight += incoming_connections[neuronid] + 1;
                }
			}*/
			free(incoming_connections);            
		}
	}
private:
    std::vector<Population*>* populations;
};

#endif
