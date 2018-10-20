#ifndef PCX_COMM_H
#define PCX_COMM_H

// Communication pattern based on Hoefler, T., Siebert, C., & Lumsdaine, A. (2010). Scalable communication protocols for dynamic sparse data exchange. ACM SIGPLAN Notices, 45(5), 159. https://doi.org/10.1145/1837853.1693476

#include <algorithm>
#include <vector>
#include "P2PCommunicator.h"
#include "Utils.h"


class PCXCommunicator : public P2PCommunicator {
public:
	
	PCXCommunicator(int p_id, int np, int master, int population_size) : P2PCommunicator(p_id,np,master,population_size) {
		recvcounts = (int*)malloc(sizeof(int) * np);
		for(int ii=0; ii < np; ii++) recvcounts[ii] = 1;
		local_table.resize(np);
		requests = (MPI_Request*) malloc(sizeof(MPI_Request) * np);
	}
	virtual ~PCXCommunicator() {
		free(recvcounts);
		free(requests);
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		comm_step++;

		// prepare list I of destinations of data
		// allocate local table with P entries, init to 0
		// for each process I need to send message, set row target(i) to 1
		// MPI_Reduce_scatter to spread the number of messages to be sent to each process
		// nonblocking send to all destinations
		// for each round (as many rounds as messages I should receive)
		// 		blocking probe for incoming message
		//		receive message
		
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		generate_2D_structure(propagate, &to_send, num_processes - 1);
		
		for(int ii=0; ii < num_processes; ii++) {
			local_table[ii] = to_send[ii].size() > 0 ? 1 : 0;
		}
		communication_sent += sizeof(int) * 2;
		num_messages+= 2; // Reduce scatter is assumed to be implemented with 2 messages per process
		int receive_count;
		
		MPI_Reduce_scatter(&local_table[0],&receive_count,recvcounts,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		
		// nonblocking synchronous send to each process that requires info
		int send_size_bytes, send_size;
		bool active_comm_step = false;
		for(int ii=0; ii < num_processes; ii++) {
			send_size = to_send[ii].size();
			if(send_size == 0) {
				requests[ii] = MPI_REQUEST_NULL;
				continue;
			}
			active_comm_step = true;
			send_size_bytes = max(sizeof(unsigned int) * send_size,0);
			/*runtime_neighbours++;
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size_bytes);
#endif
			communication_sent += send_size_bytes;
			num_messages++;	*/	
			stat_send_message_to(send_size_bytes,ii);
				
			MPI_Isend(&(to_send[ii][0]),send_size,MPI_UNSIGNED,ii,comm_step,MPI_COMM_WORLD,&requests[ii]);
			
		}
		if(active_comm_step) non_empty_comm_step++;
		
		int flag;
		MPI_Status st;
		while(receive_count > 0) {
			MPI_Iprobe(MPI_ANY_SOURCE,comm_step,MPI_COMM_WORLD,&flag,&st);
			if(flag) {
				int count;
				MPI_Get_count(&st,MPI_INT,&count);
				unsigned int* buffer = (unsigned int*)malloc(sizeof(unsigned int) * count);
				MPI_Recv(buffer,count,MPI_UNSIGNED,st.MPI_SOURCE,comm_step,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for(int ss=0; ss < count; ss++) {
					received_propagate->push_back(buffer[ss]);
				}
				free(buffer);
				receive_count--;
			}
		}
		MPI_Waitall(num_processes,requests,MPI_STATUS_IGNORE); // necessary if send messages are large (one process may have received all messages but is still sending, as it is async)
	}
	
private:
	std::vector<int> local_table;
	int* recvcounts;
	MPI_Request* requests;
};


/*class PCXCommunicator : public Communicator {
public:
	
	PCXCommunicator(int p_id, int np, int master) : Communicator(p_id,np,master) {
		
	}
	virtual ~PCXCommunicator() {}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, idx_t* partitioning, std::vector<std::vector<int> >* connections, std::vector<unsigned int>* received_propagate) {
		comm_step++;
		// prepare list I of destinations of data
		// allocate local table with P entries, init to 0
		// for each process I need to send message, set row target(i) to 1
		// MPI_Reduce_scatter to spread the number of messages to be sent to each process
		// nonblocking send to all destinations
		// for each round (as many rounds as messages I should receive)
		// 		blocking probe for incoming message
		//		receive message
		
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
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
			for(int ss=0; ss < size; ss++) {
				if(partitioning[cons[ss]] != process_id) {
					//check if the presynaptic neuron was already in the to be sent list
					if(std::find(to_send[partitioning[cons[ss]]].begin(),to_send[partitioning[cons[ss]]].end(),message) == to_send[partitioning[cons[ss]]].end()) {
						to_send[partitioning[cons[ss]]].push_back(message);
						destinations++;
					}
					if(destinations >= num_processes - 1) { // already notifying all partitions of the spike of neuron propagate[ii]
						break;
					}
				}
			}
			spikes_sent += destinations; 
		}
		
		std::vector<int> local_table(num_processes,0);
		int* recvcounts = (int*)malloc(sizeof(int) * num_processes);
		for(int ii=0; ii < num_processes; ii++) {
			if(to_send[ii].size() > 0) local_table[ii] = 1;
			recvcounts[ii] = 1;
		}
		int receive_count;
		
		MPI_Reduce_scatter(&local_table[0],&receive_count,recvcounts,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		
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
		
	}
	
	
};*/

#endif
