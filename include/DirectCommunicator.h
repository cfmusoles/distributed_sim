#ifndef DIRECT_COMM_H
#define DIRECT_COMM_H

#include <algorithm>
#include "P2PCommunicator.h"
#include "Utils.h"


class DirectCommunicator : public P2PCommunicator {
public:
	DirectCommunicator(int p_id, int np, int master, int population_size) : P2PCommunicator(p_id,np,master,population_size) {
		
	}
	virtual ~DirectCommunicator() {
		
	}
	
	virtual void send_receive(std::vector<unsigned int>* propagate, std::vector<unsigned int>* received_propagate) {
		comm_step++;
		// create 2D array of processes-spikes
		std::vector<std::vector<unsigned int> > to_send(num_processes);
		generate_2D_structure(propagate, &to_send, num_processes - 1);
		
		// send receive data
		std::vector<int> processes_to_listen;
		if(process_id != master_node) {
			// communicate target of P2P to MASTER_NODE
			std::vector<int> target_nodes;
			for(int ii=0; ii < num_processes; ii++) {
				if(to_send[ii].size() > 0) {
					target_nodes.push_back(ii);
				}
			}
			
			communication_sent += max(sizeof(int) * target_nodes.size(),0);
			MPI_Send(&target_nodes[0],target_nodes.size(),MPI_INT,master_node,0,MPI_COMM_WORLD);
			num_messages++;
			// receive how many messages to expect
			MPI_Status status;
			MPI_Probe(master_node,1,MPI_COMM_WORLD,&status);
			int count;
			MPI_Get_count(&status,MPI_INT,&count);
			int* buffer = (int*)malloc(sizeof(int) * count);
			MPI_Recv(buffer,count,MPI_INT,master_node,1,MPI_COMM_WORLD,&status);
			for(int tt=0; tt < count; tt++) {
				processes_to_listen.push_back(buffer[tt]);
			}
			free(buffer);
		} else {
			int* flags = (int*)malloc(sizeof(int) * num_processes);
			MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status) * num_processes);
			int receivs = 0;
			std::vector<std::vector<int> > target_messages(num_processes);
			while(receivs < num_processes - 1) {
				for(int ii=0; ii < num_processes; ii++) {
					if(ii == master_node) continue;
					MPI_Iprobe(ii,0,MPI_COMM_WORLD,&flags[ii],&status[ii]);
					if(flags[ii]) {
						// message received, process it
						int count;
						MPI_Get_count(&status[ii],MPI_INT,&count);
						int* buffer = (int*)malloc(sizeof(int) * count);
						MPI_Recv(buffer,count,MPI_INT,ii,0,MPI_COMM_WORLD,&status[ii]);
						for(int tt=0; tt < count; tt++) {
							// if intended target is MASTER_NODE, store locally
							if(buffer[tt] == master_node)
								processes_to_listen.push_back(ii);
							else
								target_messages[buffer[tt]].push_back(ii);
						}
						free(buffer);
						receivs++;
					}
				}
			}
			
			free(flags);
			free(status);
			// notify all processes of how many messages to expect
			// also add if they will hear from MASTER_NODE
			for(int ii=0; ii < num_processes; ii++) {
				if(ii == master_node) continue;
				if(to_send[ii].size() > 0)
					target_messages[ii].push_back(master_node);
			}
			for(int ii=0; ii < num_processes; ii++) {
				if(ii == master_node) continue;
				if(target_messages[ii].size() == 0) {
					int dummy;
					MPI_Send(&dummy,0,MPI_INT,ii,1,MPI_COMM_WORLD);
				} else {
					int send_size = max(sizeof(int) * target_messages[ii].size(),0);
					communication_sent += send_size;
					MPI_Send(&target_messages[ii][0],target_messages[ii].size(),MPI_INT,ii,1,MPI_COMM_WORLD);
				}
				num_messages++;
			}
		}
		// FOr all, master and workers
		// once all processes know which processes to listen, start P2P listening and sending
		// send
		MPI_Request* reqs = (MPI_Request*)malloc(sizeof(MPI_Request) * num_processes);
		bool active_comm_step = false;
		for(int ii=0; ii < num_processes; ii++) {
			if(to_send[ii].size() > 0) {
				//printf("%i: sends %li to %i\n",process_id,propagate[ii].size(),ii);
				active_comm_step=true;
				int send_size = max(sizeof(int) * to_send[ii].size(),0);
				/*
#ifdef ADVANCED_COMM_STATS
				messages_sent.push_back(send_size);
#endif
				runtime_neighbours++;
				communication_sent += send_size;
				num_messages++;*/
				stat_send_message_to(send_size,ii);
				MPI_Isend(&to_send[ii][0],to_send[ii].size(),MPI_INT,ii,2,MPI_COMM_WORLD,&reqs[ii]);
				MPI_Request_free(&reqs[ii]);
			}				
		}
		if(active_comm_step) non_empty_comm_step++;

		free(reqs);
		// receive
		int* flags = (int*)malloc(sizeof(int) * processes_to_listen.size());
		MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status) * processes_to_listen.size());
		int receivs = 0;
		while(receivs < processes_to_listen.size()) {
			for(int ii=0; ii < processes_to_listen.size(); ii++) {
				MPI_Iprobe(processes_to_listen[ii],2,MPI_COMM_WORLD,&flags[ii],&status[ii]);
				if(flags[ii]) {
					// process spikes
					int count;
					MPI_Get_count(&status[ii],MPI_INT,&count);
					unsigned int* buffer = (unsigned int*)malloc(sizeof(unsigned int) * count);
					MPI_Recv(buffer,count,MPI_INT,processes_to_listen[ii],2,MPI_COMM_WORLD,&status[ii]);
					// buffer contains presynaptic spiking neurons from process process_to_listen[ii]
					//printf("%i: receives %i from %i\n",process_id,count,processes_to_listen[ii]);
					for(int ss=0; ss < count; ss++) {
						// dont think checking if spike has already been received is necessary
						//if(std::find(received_propagate->begin(),received_propagate->end(),buffer[ss]) == received_propagate->end())
						received_propagate->push_back(buffer[ss]);
					}
					free(buffer);
					receivs++;
				}
			}
		}
		free(flags);
		free(status);

	}
private:
	
		
};


/*class DirectCommunicator : public Communicator {
public:
	DirectCommunicator(int p_id, int np, int master) : Communicator(p_id,np,master) {
		
	}
	virtual ~DirectCommunicator() {
		
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
		
		// send receive data
		std::vector<int> processes_to_listen;
		if(process_id != master_node) {
			// communicate target of P2P to MASTER_NODE
			std::vector<int> target_nodes;
			for(int ii=0; ii < num_processes; ii++) {
				if(to_send[ii].size() > 0)
					target_nodes.push_back(ii);
			}
			communication_sent += max(sizeof(int) * target_nodes.size(),0);
			MPI_Send(&target_nodes[0],target_nodes.size(),MPI_INT,master_node,0,MPI_COMM_WORLD);
			num_messages++;
			// receive how many messages to expect
			MPI_Status status;
			MPI_Probe(master_node,1,MPI_COMM_WORLD,&status);
			int count;
			MPI_Get_count(&status,MPI_INT,&count);
			int* buffer = (int*)malloc(sizeof(int) * count);
			MPI_Recv(buffer,count,MPI_INT,master_node,1,MPI_COMM_WORLD,&status);
			for(int tt=0; tt < count; tt++) {
				processes_to_listen.push_back(buffer[tt]);
			}
			free(buffer);
		} else {
			int* flags = (int*)malloc(sizeof(int) * num_processes);
			MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status) * num_processes);
			int receivs = 0;
			std::vector<std::vector<int> > target_messages(num_processes);
			while(receivs < num_processes - 1) {
				for(int ii=0; ii < num_processes; ii++) {
					if(ii == master_node) continue;
					MPI_Iprobe(ii,0,MPI_COMM_WORLD,&flags[ii],&status[ii]);
					if(flags[ii]) {
						// message received, process it
						int count;
						MPI_Get_count(&status[ii],MPI_INT,&count);
						int* buffer = (int*)malloc(sizeof(int) * count);
						MPI_Recv(buffer,count,MPI_INT,ii,0,MPI_COMM_WORLD,&status[ii]);
						for(int tt=0; tt < count; tt++) {
							// if intended target is MASTER_NODE, store locally
							if(buffer[tt] == master_node)
								processes_to_listen.push_back(ii);
							else
								target_messages[buffer[tt]].push_back(ii);
						}
						free(buffer);
						receivs++;
					}
				}
			}
			
			free(flags);
			free(status);
			// notify all processes of how many messages to expect
			// also add if they will hear from MASTER_NODE
			for(int ii=0; ii < num_processes; ii++) {
				if(ii == master_node) continue;
				if(to_send[ii].size() > 0)
					target_messages[ii].push_back(master_node);
			}
			for(int ii=0; ii < num_processes; ii++) {
				if(ii == master_node) continue;
				if(target_messages[ii].size() == 0) {
					int dummy;
#ifdef ADVANCED_COMM_STATS
					messages_sent.push_back(0);
#endif
					//communication_sent += 1;
					MPI_Send(&dummy,0,MPI_INT,ii,1,MPI_COMM_WORLD);
				} else {
					int send_size = max(sizeof(int) * target_messages[ii].size(),0);
#ifdef ADVANCED_COMM_STATS
					messages_sent.push_back(send_size);
#endif
					communication_sent += send_size;
					MPI_Send(&target_messages[ii][0],target_messages[ii].size(),MPI_INT,ii,1,MPI_COMM_WORLD);
				}
				num_messages++;
			}
		}
		// FOr all, master and workers
		// once all processes know which processes to listen, start P2P listening and sending
		// send
		MPI_Request* reqs = (MPI_Request*)malloc(sizeof(MPI_Request) * num_processes);
		for(int ii=0; ii < num_processes; ii++) {
			if(to_send[ii].size() > 0) {
				//printf("%i: sends %li to %i\n",process_id,propagate[ii].size(),ii);
				int send_size = max(sizeof(int) * to_send[ii].size(),0);
#ifdef ADVANCED_COMM_STATS
				messages_sent.push_back(send_size);
#endif
				communication_sent += send_size;
				MPI_Isend(&to_send[ii][0],to_send[ii].size(),MPI_INT,ii,2,MPI_COMM_WORLD,&reqs[ii]);
				num_messages++;
				MPI_Request_free(&reqs[ii]);
			}				
		}
		free(reqs);
		// receive
		int* flags = (int*)malloc(sizeof(int) * processes_to_listen.size());
		MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status) * processes_to_listen.size());
		int receivs = 0;
		while(receivs < processes_to_listen.size()) {
			for(int ii=0; ii < processes_to_listen.size(); ii++) {
				MPI_Iprobe(processes_to_listen[ii],2,MPI_COMM_WORLD,&flags[ii],&status[ii]);
				if(flags[ii]) {
					// process spikes
					int count;
					MPI_Get_count(&status[ii],MPI_INT,&count);
					unsigned int* buffer = (unsigned int*)malloc(sizeof(unsigned int) * count);
					MPI_Recv(buffer,count,MPI_INT,processes_to_listen[ii],2,MPI_COMM_WORLD,&status[ii]);
					// buffer contains presynaptic spiking neurons from process process_to_listen[ii]
					//printf("%i: receives %i from %i\n",process_id,count,processes_to_listen[ii]);
					for(int ss=0; ss < count; ss++) {
						// dont think checking if spike has already been received is necessary
						//if(std::find(received_propagate->begin(),received_propagate->end(),buffer[ss]) == received_propagate->end())
						received_propagate->push_back(buffer[ss]);
					}
					free(buffer);
					receivs++;
				}
			}
		}
		free(flags);
		free(status);

	}
private:
	
		
};*/

#endif
