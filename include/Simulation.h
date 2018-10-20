#ifndef SIMULATION__H
#define SIMULATION__H

#include <string>
#include "NeuronParams.h"
#include <vector>
#include <metis.h>

enum Population_type {
	EXC,
	INH
};

class Population {
public:
	std::string name;
	int num_neurons;
	Population_type type;
	NeuronParams* cell_params;
	std::vector<int> neuron_ids;
	
	Population(std::string nm, int n, Population_type t, NeuronParams* c) {
		num_neurons = n;
		type = t;
		name = nm;
		cell_params = c;
	}
	~Population() {}
	
	bool is_in_population(int id) {
		for(int ii=0; ii < neuron_ids.size(); ii++) {
			if(id == neuron_ids[ii]) return true;
		}
		return false;
	}
};

struct Model {
	idx_t population_size;
	int total_connections;
	int** interconnections;
	bool null_compute;
	int* interconnections_size;
	char* hypergraph_file = NULL;

	bool store_in_file;
};



#endif

