#ifndef ALL_GATHER_COMM_H
#define ALL_GATHER_COMM_H

#include <algorithm>
#include "Communicator.h"
#include "Utils.h"

class AllGatherCommunicator : public Communicator {
public:
	AllGatherCommunicator(int p_id, int np, int master, int population_size) : Communicator(p_id,np,master,population_size) {
		recvcounts = (int*)malloc(sizeof(int) * np);
		displs = (int*)malloc(sizeof(int) * np);
	}
	virtual ~AllGatherCommunicator() {
		free(recvcounts);
		free(displs);
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		comm_step++;
		int num_sends = propagate->size();
		spikes_sent += num_sends * (num_processes-1);
		if(num_sends > 0) non_empty_comm_step++;
		int send_size = sizeof(int) * 1 * (num_processes - 1);
		communication_sent += send_size;
		num_messages += (num_processes - 1);
		MPI_Allgather(&num_sends,1,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);
		int counter = 0;
		for(int ii=0; ii < num_processes; ii++) {
			displs[ii] = counter;
			counter += recvcounts[ii]; 
		}
		unsigned int* recvbuffer = (unsigned int*)malloc(sizeof(unsigned int) * counter);
		//send_size = max(sizeof(int) * propagate->size(),0) * (num_processes - 1);
		send_size = max(sizeof(int) * propagate->size(),0);
		if(send_size > 0) {
			for(int ii=0; ii < num_processes; ii++) {
				if(ii == process_id) continue;
				/*
	#ifdef ADVANCED_COMM_STATS
				messages_sent.push_back(send_size);
	#endif
				communication_sent += send_size;
				num_messages++;
				runtime_neighbours++;*/
				stat_send_message_to(send_size,ii);
			}
		}

		MPI_Allgatherv(&(*propagate)[0],propagate->size(),MPI_INT,recvbuffer,recvcounts,displs,MPI_INT,MPI_COMM_WORLD);
		for(int ii=0; ii < num_processes; ii++)  {
			if(process_id == ii)
				continue;
			for(int jj=0; jj < recvcounts[ii]; jj++) {
				unsigned int preneuron = recvbuffer[jj + displs[ii]];
				// dont think checking if spike has already been received is necessary
				//if(std::find(received_propagate->begin(),received_propagate->end(),preneuron) == received_propagate->end())
				received_propagate->push_back(preneuron);
				
			}
		}
		free(recvbuffer);
	}
private:
	int* recvcounts;
	int* displs;
		
};

#endif
