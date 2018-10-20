#ifndef SUBSCRIBER_COMM_H
#define SUBSCRIBER_COMM_H

#include <algorithm>
#include <vector>
#include "P2PCommunicator.h"
#include "Utils.h"


class SubscriberCommunicator : public P2PCommunicator {
public:
	
	SubscriberCommunicator(int p_id, int np, int master, int population_size) : P2PCommunicator(p_id,np,master,population_size) {
		
	}
	
	virtual ~SubscriberCommunicator() {
		free(requests);
	}
	
	virtual void from_to_connection(int from, int to, idx_t* partitioning) {
		P2PCommunicator::from_to_connection(from, to, partitioning);
		
		if(partitioning[from] == process_id && partitioning[to] != process_id) {
			if(std::find(data_to.begin(),data_to.end(),partitioning[to]) == data_to.end()) 
				data_to.push_back(partitioning[to]);
		}
		if(partitioning[from] != process_id && partitioning[to] == process_id) {
			if(std::find(data_from.begin(),data_from.end(),partitioning[from]) == data_from.end()) 
				data_from.push_back(partitioning[from]);
		}
	}
	
	virtual void final_setup() {
		requests = (MPI_Request*) malloc(sizeof(MPI_Request) * data_to.size());
	}
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		
		comm_step++;
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		generate_2D_structure(propagate, &to_send, data_to.size());
		
		// start P2P listening and sending
		// send
		int send_size;
		for(int ii=0; ii < data_to.size(); ii++) {
			int dest = data_to[ii];
			send_size = max(sizeof(int) * to_send[dest].size(),0);
			/*
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size);
#endif
			communication_sent += send_size;
			num_messages++;
			runtime_neighbours++;*/
			stat_send_message_to(send_size,dest);
			MPI_Isend(&to_send[dest][0],to_send[dest].size(),MPI_INT,dest,comm_step,MPI_COMM_WORLD,&requests[ii]);
		}
		non_empty_comm_step++; // always sends messages to neighbours

		// receive
		int receivs = 0;
		MPI_Status st;
		while(receivs < data_from.size()) {
			MPI_Probe(MPI_ANY_SOURCE,comm_step,MPI_COMM_WORLD,&st);
			// process spikes
			int count;
			MPI_Get_count(&st,MPI_INT,&count);
			unsigned int* buffer = (unsigned int*)malloc(sizeof(unsigned int) * count);
			MPI_Recv(buffer,count,MPI_INT,st.MPI_SOURCE,comm_step,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			// buffer contains presynaptic spiking neurons from process st.MPI_SOURCE
			//printf("%i: receives %i from %i\n",process_id,count,st.MPI_SOURCE);
			for(int ss=0; ss < count; ss++) {
				received_propagate->push_back(buffer[ss]);
			}
			free(buffer);
			receivs++;
		}	
		
		MPI_Waitall(data_to.size(),requests,MPI_STATUS_IGNORE); // necessary if send messages are large (one process may have received all messages but is still sending, as it is async)
	}
	
	
private:
	// each process should know which other processes it needs to communicate data to and from
	std::vector<int> data_from;
	std::vector<int> data_to;
	MPI_Request* requests;
};



/*class SubscriberCommunicator : public Communicator {
public:
	
	SubscriberCommunicator(int p_id, int np, int master, int population_size) : Communicator(p_id,np,master) {
		connection_to.resize(population_size);
	}
	
	virtual ~SubscriberCommunicator() {}
	
	virtual void from_to_connection(int from, int to, idx_t* partitioning) {
		if(partitioning[from] == process_id && partitioning[to] != process_id) {
			if(std::find(data_to.begin(),data_to.end(),partitioning[to]) == data_to.end()) 
				data_to.push_back(partitioning[to]);
		}
		if(partitioning[from] != process_id && partitioning[to] == process_id) {
			if(std::find(data_from.begin(),data_from.end(),partitioning[from]) == data_from.end()) 
				data_from.push_back(partitioning[from]);
		}
		if(partitioning[from] == process_id && partitioning[to] != process_id) {
			if(std::find(connection_to[from].begin(),connection_to[from].end(),partitioning[to]) == connection_to[from].end()) 
				connection_to[from].push_back(partitioning[to]);
		}
	}
	
	virtual int get_sync_group_size() {
		return data_to.size() + data_from.size();
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, idx_t* partitioning, std::vector<std::vector<int> >* connections, std::vector<unsigned int>* received_propagate) {
		
		comm_step++;
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		for(int ii=0; ii < propagate->size(); ii++) {
			int neuronid, timing, part_dest, cons_size;
			unsigned int message = (*propagate)[ii];
			int destinations = 0;
			decodeSpike(message,&timing,&neuronid);
			bool current_destinations[num_processes];
			memset(current_destinations,0,num_processes*sizeof(bool));
			cons_size = connection_to[neuronid].size();
			for(int jj=0; jj < cons_size; jj++) {
				part_dest = connection_to[neuronid][jj];
				if(!current_destinations[part_dest]) {
					to_send[part_dest].push_back(message);
					current_destinations[part_dest] = true;
					destinations++;
				}
				if(destinations >= data_to.size()) { // already notifying all partitions of the spike of neuron propagate[ii]
					break;
				}
			}
			spikes_sent += destinations; 
		}
		
		// start P2P listening and sending
		// send
		int send_size;
		MPI_Request r;
		for(int ii=0; ii < data_to.size(); ii++) {
			int dest = data_to[ii];
			send_size = max(sizeof(int) * to_send[dest].size(),0);
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size);
#endif
			communication_sent += send_size;
			MPI_Isend(&to_send[dest][0],to_send[dest].size(),MPI_INT,dest,comm_step,MPI_COMM_WORLD,&r);
			num_messages++;
		}
		
		// receive
		int receivs = 0;
		MPI_Status st;
		while(receivs < data_from.size()) {
			MPI_Probe(MPI_ANY_SOURCE,comm_step,MPI_COMM_WORLD,&st);
			// process spikes
			int count;
			MPI_Get_count(&st,MPI_INT,&count);
			unsigned int* buffer = (unsigned int*)malloc(sizeof(unsigned int) * count);
			MPI_Recv(buffer,count,MPI_INT,st.MPI_SOURCE,comm_step,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			// buffer contains presynaptic spiking neurons from process st.MPI_SOURCE
			//printf("%i: receives %i from %i\n",process_id,count,st.MPI_SOURCE);
			for(int ss=0; ss < count; ss++) {
				received_propagate->push_back(buffer[ss]);
			}
			free(buffer);
			receivs++;
		}	
		
	}
	
	
private:
	// each process should know which other processes it needs to communicate data to and from
	std::vector<int> data_from;
	std::vector<int> data_to;
	
	std::vector<std::vector<int> > connection_to;
};*/


#endif
