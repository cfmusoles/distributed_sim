#ifndef RSX_COMM_H
#define RSX_COMM_H

// Communication pattern based on Hoefler, T., Siebert, C., & Lumsdaine, A. (2010). Scalable communication protocols for dynamic sparse data exchange. ACM SIGPLAN Notices, 45(5), 159. https://doi.org/10.1145/1837853.1693476

#include <algorithm>
#include <vector>
#include "P2PCommunicator.h"
#include "Utils.h"

class RSXCommunicator : public P2PCommunicator {
public:
	
	RSXCommunicator(int p_id, int np, int master,int population_size) : P2PCommunicator(p_id,np,master,population_size) {
		buffer_size = 1;
		MPI_Alloc_mem(buffer_size * sizeof(unsigned int),MPI_INFO_NULL,&winbuffer);
		memset(winbuffer,0,sizeof(unsigned int) * buffer_size);	
		MPI_Win_create(winbuffer,buffer_size * sizeof(unsigned int),sizeof(unsigned int),MPI_INFO_NULL,MPI_COMM_WORLD,&winobject);
	}
	virtual ~RSXCommunicator() {
		MPI_Free_mem(winbuffer);
		MPI_Win_free(&winobject);
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		comm_step++;
		// prepare list I of destinations of data
		// Remote memory access to notify all receiving processes
		// fence to synchronise
		// post non-blocking receives
		// send data to processes
		
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		generate_2D_structure(propagate, &to_send, num_processes - 1);
		
		MPI_Win_fence(MPI_MODE_NOPRECEDE,winobject);
		int dummy = 1;
		for(int ii=0; ii < num_processes; ii++) {
			if(to_send[ii].size() > 0) {
				MPI_Accumulate(&dummy,1,MPI_INT,ii,0,1,MPI_INT,MPI_SUM,winobject);
			}
		}
		MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE,winobject);
		
		int receive_count = winbuffer[0];
		
		// nonblocking synchronous send to each process that requires info
		MPI_Request request;
		int send_size_bytes, send_size;
		bool active_comm_step = false;
		for(int ii=0; ii < num_processes; ii++) {
			send_size = to_send[ii].size();
			if(send_size == 0) continue;
			active_comm_step = true;
			send_size_bytes = max(sizeof(int) * send_size,0);
			/*
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size_bytes);
#endif
			communication_sent += send_size_bytes;
			num_messages++;
			runtime_neighbours++;*/
			stat_send_message_to(send_size_bytes,ii);
			MPI_Isend(&(to_send[ii][0]),send_size,MPI_INT,ii,comm_step,MPI_COMM_WORLD,&request);
		}
		if(active_comm_step) non_empty_comm_step++;
		
		MPI_Status st;
		for(int ii =0; ii < receive_count; ii++) {
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
		}
		
		winbuffer[0] = 0;
	}
	
private:
	MPI_Win winobject;
	unsigned int* winbuffer;
	int buffer_size;
	
	
};


/*// TODO: Shows segfault (11) error in ShARC
class RSXCommunicator : public Communicator {
public:
	
	RSXCommunicator(int p_id, int np, int master) : Communicator(p_id,np,master) {
		buffer_size = 1;
		MPI_Alloc_mem(buffer_size * sizeof(unsigned int),MPI_INFO_NULL,&winbuffer);
		memset(winbuffer,UINT_MAX,sizeof(unsigned int) * buffer_size);	
		MPI_Win_create(winbuffer,buffer_size * sizeof(unsigned int),sizeof(unsigned int),MPI_INFO_NULL,MPI_COMM_WORLD,&winobject);
	}
	virtual ~RSXCommunicator() {
		MPI_Free_mem(winbuffer);
		MPI_Win_free(&winobject);
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, idx_t* partitioning, std::vector<std::vector<int> >* connections, std::vector<unsigned int>* received_propagate) {
		comm_step++;
		// prepare list I of destinations of data
		// Remote memory access to notify all receiving processes
		// fence to synchronise
		// post non-blocking receives
		// send data to processes
		
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
			//	if(partitioning[cons[ss]] != process_id) {
			//		//check if the presynaptic neuron was already in the to be sent list
			//		if(std::find(to_send[partitioning[cons[ss]]].begin(),to_send[partitioning[cons[ss]]].end(),message) == to_send[partitioning[cons[ss]]].end()) {
			//			to_send[partitioning[cons[ss]]].push_back(message);
			//			destinations++;
			//		}
			//		if(destinations >= num_processes - 1) { // already notifying all partitions of the spike of neuron propagate[ii]
			//			break;
			//		}
			//	}
			//}
			spikes_sent += destinations; 
		}
		
		MPI_Win_fence(MPI_MODE_NOPRECEDE,winobject);
		int dummy = 1;
		for(int ii=0; ii < num_processes; ii++) {
			if(to_send[ii].size() > 0) {
				MPI_Accumulate(&dummy,1,MPI_INT,ii,0,1,MPI_INT,MPI_SUM,winobject);
			}
		}
		MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE,winobject);
		
		int receive_count = winbuffer[0];
		
		// nonblocking synchronous send to each process that requires info
		MPI_Request request;
		int send_size_bytes, send_size;
		for(int ii=0; ii < num_processes; ii++) {
			send_size = to_send[ii].size();
			if(send_size == 0) continue;
			send_size_bytes = max(sizeof(int) * send_size,0);
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size_bytes);
#endif
			communication_sent += send_size_bytes;
			num_messages++;
			MPI_Isend(&(to_send[ii][0]),send_size,MPI_INT,ii,comm_step,MPI_COMM_WORLD,&request);
		}
		
		MPI_Status st;
		for(int ii =0; ii < receive_count; ii++) {
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
		}
		
		winbuffer[0] = 0;
	}
	
private:
	MPI_Win winobject;
	unsigned int* winbuffer;
	int buffer_size;
	
	
};*/

#endif
