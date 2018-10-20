#ifndef RMA_COMMUNICATOR__H
#define RMA_COMMUNICATOR__H

#include "P2PCommunicator.h"
#include <mpi.h>
#include <stdio.h>
#include <limits.h>

// BROKEN: possible error Integer conversion resulted in a change of sign.
// BROKEN: pop_size no longer represents the size of the global population (only local neurons)
class RMACommunicator : public P2PCommunicator {
public:
	
	RMACommunicator(int p_id, int np, int master, int pop_size) : P2PCommunicator(p_id,np,master,pop_size) {
		buffer_size = pop_size;
		MPI_Alloc_mem(buffer_size * sizeof(unsigned int),MPI_INFO_NULL,&winbuffer);
		memset(winbuffer,UINT_MAX,sizeof(unsigned int) * buffer_size);	// possible error Integer conversion resulted in a change of sign.
		MPI_Win_create(winbuffer,buffer_size * sizeof(unsigned int),sizeof(unsigned int),MPI_INFO_NULL,MPI_COMM_WORLD,&winobject);
	}
	
	virtual ~RMACommunicator() {
		
		MPI_Free_mem(winbuffer);
		MPI_Win_free(&winobject);
		MPI_Group_free(&sync_group);
	}

	virtual void from_to_connection(int from, int to, idx_t* partitioning) {
		P2PCommunicator::from_to_connection(from, to, partitioning);
		
		if(partitioning[from] == process_id && partitioning[to] != process_id) {
			if(std::find(ranks.begin(),ranks.end(),partitioning[to]) == ranks.end()) 
				ranks.push_back(partitioning[to]);
		}
		if(partitioning[from] != process_id && partitioning[to] == process_id) {
			if(std::find(ranks.begin(),ranks.end(),partitioning[from]) == ranks.end()) 
				ranks.push_back(partitioning[from]);
		}
	}
	
	virtual void final_setup() {
		int nInGroup = ranks.size();
		MPI_Group group;
		MPI_Win_get_group(winobject,&group);
		MPI_Group_incl(group,nInGroup,&ranks[0],&sync_group);
		MPI_Group_free(&group);
	}
	
	virtual int get_sync_group_size() {
		return ranks.size();
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		generate_2D_structure(propagate, &to_send, num_processes - 1);
		
		// RMA operations must be wrapped within fences
		MPI_Win_post(sync_group,0,winobject);
		MPI_Win_start(sync_group,0,winobject);
		//MPI_Win_fence(MPI_MODE_NOPRECEDE,winobject);
		
		int displ = buffer_size / num_processes * process_id;
		bool active_comm_step = false;
		for(int ii=0; ii < num_processes; ii++) {
			if(ii == process_id) continue;
			int size = to_send[ii].size();
			if(size == 0) continue;
			active_comm_step = true;
			MPI_Put(&(to_send[ii][0]),size,MPI_INT,ii,displ,size,MPI_INT,winobject);
			int send_size = max(sizeof(int) * size,0);
			/*
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size);
#endif
			communication_sent += send_size;
			num_messages++;
			runtime_neighbours++;
			*/
			stat_send_message_to(send_size,ii);
		}
		if(active_comm_step) non_empty_comm_step++;
		
		MPI_Win_complete(winobject);
		MPI_Win_wait(winobject);
		//MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE,winobject);
		// read from buffer and put spikes into received_propagate vector
		// potentially slow to scale!
		int block = buffer_size / num_processes;
		for(int ii=0; ii < num_processes; ii++) {
			if(ii == process_id) continue;
			int counter = 0;
			unsigned int* pointer = &winbuffer[block * ii];
			while(counter < block && *pointer != UINT_MAX) {
				received_propagate->push_back(*pointer);
				pointer++;
				counter++;
			}
		}
			
		// clear buffer for next iteration
		memset(winbuffer,UINT_MAX,sizeof(unsigned int) * buffer_size); // possible error Integer conversion resulted in a change of sign.
		
	}
	
private:

	MPI_Win winobject;
	unsigned int* winbuffer;
	int buffer_size;
	MPI_Group sync_group;
	std::vector<int> ranks;
};


// TODO: hangs on ShARC?
/*class RMACommunicator : public Communicator {
public:
	
	RMACommunicator(int p_id, int np, int master, int pop_size) : Communicator(p_id,np,master) {
		buffer_size = pop_size;
		MPI_Alloc_mem(buffer_size * sizeof(unsigned int),MPI_INFO_NULL,&winbuffer);
		memset(winbuffer,UINT_MAX,sizeof(unsigned int) * buffer_size);	
		MPI_Win_create(winbuffer,buffer_size * sizeof(unsigned int),sizeof(unsigned int),MPI_INFO_NULL,MPI_COMM_WORLD,&winobject);
	}
	
	virtual ~RMACommunicator() {
		
		MPI_Free_mem(winbuffer);
		MPI_Win_free(&winobject);
		MPI_Group_free(&sync_group);
	}

	virtual void from_to_connection(int from, int to, idx_t* partitioning) {
		if(partitioning[from] == process_id && partitioning[to] != process_id) {
			if(std::find(ranks.begin(),ranks.end(),partitioning[to]) == ranks.end()) 
				ranks.push_back(partitioning[to]);
		}
		if(partitioning[from] != process_id && partitioning[to] == process_id) {
			if(std::find(ranks.begin(),ranks.end(),partitioning[from]) == ranks.end()) 
				ranks.push_back(partitioning[from]);
		}
	}
	
	virtual void final_setup() {
		int nInGroup = ranks.size();
		MPI_Group group;
		MPI_Win_get_group(winobject,&group);
		MPI_Group_incl(group,nInGroup,&ranks[0],&sync_group);
		MPI_Group_free(&group);
	}
	
	virtual int get_sync_group_size() {
		return ranks.size();
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, idx_t* partitioning, std::vector<std::vector<int> >* connections, std::vector<unsigned int>* received_propagate) {
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		int part_dest;
		for(int ii=0; ii < propagate->size(); ii++) {
			int neuronid, timing;
			unsigned int message = (*propagate)[ii];
			decodeSpike(message,&timing,&neuronid);
			int size = (*connections)[neuronid].size();
			int destinations = 0;
			int* cons;
			if(size > 0) {
				cons = &((*connections)[neuronid][0]);
			} 
			
			bool current_destinations[num_processes];
			memset(current_destinations,0,num_processes*sizeof(bool));
			for(int ss=0; ss < size; ss++) {
				part_dest = partitioning[cons[ss]];
				if(part_dest != process_id) {
					//check if the presynaptic neuron was already in the to be sent list
					if(!current_destinations[part_dest]) {
						to_send[part_dest].push_back(message);
						destinations++;
						current_destinations[part_dest] = true;
					}
					if(destinations >= num_processes - 1) { // already notifying all partitions of the spike of neuron propagate[ii]
						break;
					}
				}
			}
			
			//for(int ss=0; ss < size; ss++) {
			//	if(partitioning[(*connections)[neuronid][ss]] != process_id) {
			//		//check if the presynaptic neuron was already in the to be sent list
			//		if(std::find(to_send[partitioning[(*connections)[neuronid][ss]]].begin(),to_send[partitioning[(*connections)[neuronid][ss]]].end(),(*propagate)[ii]) == to_send[partitioning[(*connections)[neuronid][ss]]].end()) {
			//			to_send[partitioning[(*connections)[neuronid][ss]]].push_back((*propagate)[ii]);
			//			destinations++;
			//		}
			//		if(destinations >= num_processes - 1) { // already notifying all partitions of the spike of neuron propagate[ii]
			//			break;
			//		}
			//	}
			//}
			spikes_sent += destinations; 
		}
		
		
		// RMA operations must be wrapped within fences
		MPI_Win_post(sync_group,0,winobject);
		MPI_Win_start(sync_group,0,winobject);
		//MPI_Win_fence(MPI_MODE_NOPRECEDE,winobject);
		
		int displ = buffer_size / num_processes * process_id;
		int total_sync = 0;
		for(int ii=0; ii < num_processes; ii++) {
			if(ii == process_id) continue;
			int size = to_send[ii].size();
			if(size == 0) continue;
			total_sync++;
			MPI_Put(&(to_send[ii][0]),size,MPI_INT,ii,displ,size,MPI_INT,winobject);
			int send_size = max(sizeof(int) * size,0);
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size);
#endif
			communication_sent += send_size;
			num_messages++;
			
		}
		MPI_Win_complete(winobject);
		MPI_Win_wait(winobject);
		//MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE,winobject);
		// read from buffer and put spikes into received_propagate vector
		// potentially slow to scale!
		int block = buffer_size / num_processes;
		for(int ii=0; ii < num_processes; ii++) {
			if(ii == process_id) continue;
			int counter = 0;
			unsigned int* pointer = &winbuffer[block * ii];
			while(counter < block && *pointer != UINT_MAX) {
				received_propagate->push_back(*pointer);
				pointer++;
				counter++;
			}
		}
			
		// clear buffer for next iteration
		memset(winbuffer,UINT_MAX,sizeof(unsigned int) * buffer_size);
		
	}
	
private:

	MPI_Win winobject;
	unsigned int* winbuffer;
	int buffer_size;
	MPI_Group sync_group;
	std::vector<int> ranks;
};*/



#endif
