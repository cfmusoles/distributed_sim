#ifndef PEX_COMM_H
#define PEX_COMM_H

// Communication pattern based on Hoefler, T., Siebert, C., & Lumsdaine, A. (2010). Scalable communication protocols for dynamic sparse data exchange. ACM SIGPLAN Notices, 45(5), 159. https://doi.org/10.1145/1837853.1693476

#include <algorithm>
#include <vector>
#include "P2PCommunicator.h"
#include "Utils.h"


class PEXCommunicator : public P2PCommunicator {
public:
	
	PEXCommunicator(int p_id, int np, int master, int population_size) : P2PCommunicator(p_id,np,master,population_size) {
		datasizebuffer = (int*)calloc(np,sizeof(int));
        recvdatasizebuffer = (int*)calloc(np,sizeof(int));
        rdispls = (int*)calloc(np,sizeof(int));
        sdispls = (int*)calloc(np,sizeof(int));
	}
	virtual ~PEXCommunicator() {
		free(datasizebuffer);
        free(recvdatasizebuffer);
        free(sdispls);
        free(rdispls);
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		comm_step++;

		// prepare list I of destinations of data
		// each process sends data sizes in a vector of P elements to each process (MPI_Alltoall)
        // use MPI_Alltoall to send / receive irregular data
        // finish when all sends and receives are done
        
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		generate_2D_structure(propagate, &to_send, num_processes - 1);
		
        std::vector<unsigned int> to_send_flat;
        int total = 0;
        int send_size_bytes;
        for(int ii=0; ii < num_processes; ii++) {
            if(ii == process_id) continue;
            int send_size = to_send[ii].size();
            datasizebuffer[ii] = send_size;  
            sdispls[ii] = total;
            total += send_size;
            to_send_flat.insert(to_send_flat.end(),to_send[ii].begin(),to_send[ii].end()); 
            if(send_size == 0) continue;
            send_size_bytes = max(sizeof(int) * send_size,0);
            /*
#ifdef ADVANCED_COMM_STATS
			messages_sent.push_back(send_size_bytes);
#endif
            communication_sent += send_size_bytes;
		    num_messages++;  
            runtime_neighbours++; */
            stat_send_message_to(send_size_bytes,ii);
        }
        communication_sent += sizeof(int) * (num_processes-1);
		num_messages+= num_processes-1;
        MPI_Alltoall(datasizebuffer,1,MPI_INT,recvdatasizebuffer,1,MPI_INT,MPI_COMM_WORLD);
        if(propagate->size() > 0) non_empty_comm_step++;

        int total_recvs = 0;
        for(int ii=0; ii < num_processes; ii++) {
            rdispls[ii] = total_recvs;
            total_recvs += recvdatasizebuffer[ii];
        }
        unsigned int* recvbuffer = (unsigned int*)malloc(total_recvs * sizeof(unsigned int));
        MPI_Alltoallv(&to_send_flat[0],datasizebuffer,sdispls,MPI_UNSIGNED,recvbuffer,recvdatasizebuffer,rdispls,MPI_UNSIGNED,MPI_COMM_WORLD);
        
        for(int ii=0; ii < num_processes; ii++)  {
			if(process_id == ii)
				continue;
			for(int jj=0; jj < recvdatasizebuffer[ii]; jj++) {
				unsigned int preneuron = recvbuffer[jj + rdispls[ii]];
				received_propagate->push_back(preneuron);
				
			}
		}
        free(recvbuffer);
        
	}
	
private:
    int* datasizebuffer;
    int* recvdatasizebuffer;
    int* rdispls;
    int* sdispls;
};

#endif