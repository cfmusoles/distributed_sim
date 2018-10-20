#ifndef NBX_COMM_H
#define NBX_COMM_H

// Communication pattern based on Hoefler, T., Siebert, C., & Lumsdaine, A. (2010). Scalable communication protocols for dynamic sparse data exchange. ACM SIGPLAN Notices, 45(5), 159. https://doi.org/10.1145/1837853.1693476

#include <mpi.h>
#include <algorithm>
#include <vector>
#include "P2PCommunicator.h"
#include "Utils.h"


class NBXCommunicator : public P2PCommunicator {

public:
	
	NBXCommunicator(int p_id, int np, int master, int population_size) : P2PCommunicator(p_id,np,master,population_size) {
		requests = (MPI_Request*) malloc(sizeof(MPI_Request) * np);
		
	}
	virtual ~NBXCommunicator() {
		free(requests);
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		comm_step++;
		// prepare list I of destinations of data
		// for each i in I, send message to destination (nonblocking synchronous comm: MPI_Issend)
		
		// while not done
		// 		nonblocking probe to check for incoming messages (any source)
		// 			receive if detected
		//		if need_to_check_barrier
		//			done = ibarrier has completed
		// 		else 
		// 			need_to_check_barrier = all sends finished
		
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		generate_2D_structure(propagate, &to_send, num_processes-1);
		
		// nonblocking synchronous send to each process that requires info
		std::vector<MPI_Request*> requests_list;
		int send_size_bytes, send_size;
		bool active_comm_step = false;
		for(int ii=0; ii < num_processes; ii++) {
			send_size = to_send[ii].size();
			if(send_size == 0) continue;
			active_comm_step = true;
			send_size_bytes = max(sizeof(int) * send_size,0);
			/*runtime_neighbours++;
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size_bytes);
#endif
			communication_sent += send_size_bytes;
			num_messages++;*/
			stat_send_message_to(send_size_bytes,ii);
			MPI_Issend(&(to_send[ii][0]),send_size,MPI_UNSIGNED,ii,comm_step,MPI_COMM_WORLD,&requests[ii]);
			requests_list.push_back(&requests[ii]);
		}
		if(active_comm_step) non_empty_comm_step++;

		bool do_work = true;
		bool check_barrier = false;
		MPI_Request ibarrier_request;
		int flag;
		MPI_Status status;
		
		do {
			// check for any incoming message
			MPI_Iprobe(MPI_ANY_SOURCE,comm_step,MPI_COMM_WORLD,&flag,&status);
			if(flag) {
				// read and receive message
				int count;
				MPI_Get_count(&status,MPI_INT,&count);
				unsigned int* buffer = (unsigned int*)malloc(sizeof(unsigned int) * count);
				MPI_Recv(buffer,count,MPI_UNSIGNED,status.MPI_SOURCE,comm_step,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				// buffer contains presynaptic spiking neurons from process st.MPI_SOURCE
				//printf("%i: receives %i from %i\n",process_id,count,st.MPI_SOURCE);
				for(int ss=0; ss < count; ss++) {
					received_propagate->push_back(buffer[ss]);
				}
				free(buffer);
			}
			
			if(check_barrier) {
				// check if MPI_Ibarrier has been reached by all, if so, do_work = false
				MPI_Test(&ibarrier_request, &flag,MPI_STATUS_IGNORE);
				do_work = !flag;
			} else {
				// check if all MPI_Issend have completed, if so check_barrier = true
				int request_size = requests_list.size();
				int requests_left = request_size;
				for(int ii=0; ii < request_size; ii++) {
					if(*requests_list[ii] == MPI_REQUEST_NULL) {
						requests_left--;
						continue;
					}
					MPI_Test(requests_list[ii],&flag,MPI_STATUS_IGNORE);
					if(flag) {
						requests_left--;
					}			
				}
				if(requests_left == 0) {
					MPI_Ibarrier(MPI_COMM_WORLD,&ibarrier_request);
					check_barrier = true;
				}
			}
		} while(do_work);
		
		
	}
private:
	MPI_Request* requests;
	
};


/*class NBXCommunicator : public Communicator {
public:
	
	NBXCommunicator(int p_id, int np, int master) : Communicator(p_id,np,master) {
		
	}
	virtual ~NBXCommunicator() {}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, idx_t* partitioning, std::vector<std::vector<int> >* connections, std::vector<unsigned int>* received_propagate) {
		comm_step++;
		// prepare list I of destinations of data
		// for each i in I, send message to destination (nonblocking synchronous comm: MPI_Issend)
		
		// while not done
		// 		nonblocking probe to check for incoming messages (any source)
		// 			receive if detected
		//		if need_to_check_barrier
		//			done = ibarrier has completed
		// 		else 
		// 			need_to_check_barrier = all sends finished
		
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
		
		// nonblocking synchronous send to each process that requires info
		MPI_Request* requests = (MPI_Request*) malloc(sizeof(MPI_Request) * num_processes);
		std::vector<MPI_Request*> requests_list;
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
			MPI_Issend(&(to_send[ii][0]),send_size,MPI_INT,ii,comm_step,MPI_COMM_WORLD,&requests[ii]);
			requests_list.push_back(&requests[ii]);
		}
		
		bool do_work = true;
		bool check_barrier = false;
		MPI_Request ibarrier_request;
		int flag;
		MPI_Status status;
		
		do {
			// check for any incoming message
			MPI_Iprobe(MPI_ANY_SOURCE,comm_step,MPI_COMM_WORLD,&flag,&status);
			if(flag) {
				// read and receive message
				int count;
				MPI_Get_count(&status,MPI_INT,&count);
				unsigned int* buffer = (unsigned int*)malloc(sizeof(unsigned int) * count);
				MPI_Recv(buffer,count,MPI_INT,status.MPI_SOURCE,comm_step,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				// buffer contains presynaptic spiking neurons from process st.MPI_SOURCE
				//printf("%i: receives %i from %i\n",process_id,count,st.MPI_SOURCE);
				for(int ss=0; ss < count; ss++) {
					received_propagate->push_back(buffer[ss]);
				}
				free(buffer);
			}
			
			if(check_barrier) {
				// check if MPI_Ibarrier has been reached by all, if so, do_work = false
				MPI_Test(&ibarrier_request, &flag,MPI_STATUS_IGNORE);
				do_work = !flag;
			} else {
				// check if all MPI_Issend have completed, if so check_barrier = true
				int request_size = requests_list.size();
				int requests_left = request_size;
				for(int ii=0; ii < request_size; ii++) {
					if(*requests_list[ii] == MPI_REQUEST_NULL) {
						requests_left--;
						continue;
					}
					MPI_Test(requests_list[ii],&flag,MPI_STATUS_IGNORE);
					if(flag) {
						requests_left--;
					}			
				}
				if(requests_left == 0) {
					MPI_Ibarrier(MPI_COMM_WORLD,&ibarrier_request);
					check_barrier = true;
				}
			}
		} while(do_work);
		
		free(requests);
		
	}
	
	
};*/

#endif
