#ifndef RANDOM_CONNECTIVITY__H
#define RANDOM_CONNECTIVITY__H

#include "Connectivity.h"
#include "Utils.h"
#include <random>
#include <iomanip>
#include <map>
#include <stdio.h>

class RandomConnectivity : public Connectivity {
public:
	float p_conn;
	
	RandomConnectivity(Population* f, Population* t, float p) : Connectivity(f,t) {
		p_conn = p;
	}
	virtual ~RandomConnectivity() {}
	
	virtual void generate_connectivity(int pop_size) {
		std::default_random_engine generator;
		
		generator.seed(rand01());
		std::binomial_distribution<> d(to->num_neurons, p_conn);

		for(int ff=0; ff < from->num_neurons; ff++) {
			int num_connections = d(generator);
			int rand_conn = rand01() * (float)to->num_neurons;
			
			// NEW VERSION (2.25x faster optimisation): each population selects random starting target
			// in 'to' population, and connects to the subsequent neurons (which are randomised in order)
			if(num_connections == 0) continue;
			connections[ff].resize(num_connections);
			if(num_connections + rand_conn > to->num_neurons) {
				int diff = to->num_neurons - rand_conn;
				memcpy(&connections[ff][0],&to->neuron_ids[rand_conn],sizeof(int) * diff);
				num_connections = num_connections-diff;
				memcpy(&connections[ff][diff],&to->neuron_ids[0],sizeof(int) * num_connections);
			} else {
				memcpy(&connections[ff][0],&to->neuron_ids[rand_conn],sizeof(int) * num_connections);
			}
			
			
			// OLD VERSION: each neuron randomly selects target neurons in 'to' population
			/*connections[ff].reserve(num_connections);
			for(int tt =0; tt < num_connections; tt++) {
				int rand_conn = 0;
				do {
					rand_conn = rand01() * (float)to->num_neurons;
				} while(std::find(connections[ff].begin(), connections[ff].end(), to->neuron_ids[rand_conn]) != connections[ff].end());
				connections[ff].push_back(to->neuron_ids[rand_conn]);
			}*/
		}
	}
	
};

#endif
