#ifndef P2P_COMM_H
#define P2P_COMM_H

#include <algorithm>
#include "Communicator.h"
#include "Utils.h"

class P2PCommunicator : public Communicator {
protected:
	//std::vector<std::vector<int> > neuron_connection_to;
	
	void generate_2D_structure(std::vector<unsigned int>* propagate, std::vector<std::vector<unsigned int> >* to_send, int send_limit) {
		int ii, prop_size = propagate->size();
        if(prop_size == 0) return;
		unsigned int* prop_value;
		for(ii=0, prop_value = &propagate->at(0); ii < prop_size; prop_value++, ii++) {
			int neuronid, timing, cons_size;
			unsigned int message = *prop_value;
			decodeSpike(message,&timing,&neuronid);
			cons_size = neuron_connection_to[neuronid].size();
			int* part_dest;
			int jj;
			for(part_dest = &neuron_connection_to[neuronid][0], jj=0; jj < cons_size; part_dest++, jj++) {
				to_send->at(*part_dest).push_back(message);
			}
			spikes_sent += cons_size; 
			
		}
	}

public:
	P2PCommunicator(int p_id, int np, int master, int population_size) : Communicator(p_id,np,master,population_size) {
		//neuron_connection_to.resize(population_size);
	}
	virtual ~P2PCommunicator() { }
	
	/*virtual void from_to_connection(int from, int to, idx_t* partitioning) {
		Communicator::from_to_connection(from,to,partitioning);
		if(partitioning[from] == process_id && partitioning[to] != process_id) {
			if(std::find(neuron_connection_to[from].begin(),neuron_connection_to[from].end(),partitioning[to]) == neuron_connection_to[from].end()) 
				neuron_connection_to[from].push_back(partitioning[to]);
		}
	}*/
	
		
};

#endif
