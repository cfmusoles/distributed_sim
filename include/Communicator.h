#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include "Simulation.h"
#include <metis.h>
#include <set>
#include <algorithm>

// Interface for objects that want to become Communicators between MPI processes
// They are responsible for sharing data (spike info) between processes
class Communicator {
public:
	//stats
	long int communication_sent;
	long int num_messages;
	long int spikes_sent;
	std::vector<std::vector<int> > messages_sent;
	int comm_step;
	int non_empty_comm_step;
	std::set<int> physical_neighbours;
	long int runtime_neighbours;
	std::vector<std::vector<int> > neuron_connection_to;
	

	Communicator(int p_id, int np, int master, int population_size) {
		process_id = p_id;
		num_processes = np;
		master_node = master;
		//stats
		communication_sent = 0;
		num_messages = 0;
		spikes_sent = 0;
		comm_step = 0;
		runtime_neighbours = 0;
		non_empty_comm_step = 0;
		neuron_connection_to.resize(population_size);
		messages_sent.resize(np);
	}
	virtual ~Communicator() { }

	// notification to communicator of a synaptic connection from -> to neurons
	virtual void from_to_connection(int from, int to, idx_t* partitioning) {
		if(partitioning[from] == process_id && partitioning[to] != process_id)	physical_neighbours.insert(partitioning[to]);
		if(partitioning[to] == process_id && partitioning[from] != process_id)	physical_neighbours.insert(partitioning[from]);

		// store list of partitions that have a connection coming from each pre-synaptic neuron
		if(partitioning[from] == process_id && partitioning[to] != process_id) {
			if(std::find(neuron_connection_to[from].begin(),neuron_connection_to[from].end(),partitioning[to]) == neuron_connection_to[from].end()) {
				neuron_connection_to[from].push_back(partitioning[to]);
			}
		}
	}
	// final setup operations before start simulation
	virtual void final_setup() { }
	// perform send - receive across MPI processes
	//virtual void send_receive(std::vector<unsigned int>* propagate, idx_t* partitioning, std::vector<std::vector<int> >* connections, std::vector<unsigned int>* received_propagate) {}
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {}
	// get synchronisation group size (number of processes to send / receive info from / to)
	virtual int get_sync_group_size() {
		return physical_neighbours.size();
	}
	
protected:
	int process_id;
	int num_processes;
	int master_node;

	// for statistics
	// Only for data transfers to destination processes (spikes, not comm overhead)
	void stat_send_message_to(int size, int destination) {
#if defined(ADVANCED_COMM_STATS) || defined(ADVANCED_COMM_STATS_MATRIX_ONLY)
		messages_sent[destination].push_back(size);
#endif
		communication_sent += size;
		num_messages++;
		runtime_neighbours++;
	}
	
};

#endif
