#ifndef PARALLEL_PARTITION_PARTITIONING
#define PARALLEL_PARTITION_PARTITIONING

#include "parmetis.h"
#include <vector>
#include <set>
#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"

class ParallelPartitionPartitioning : public Partitioning {
public:
			
	ParallelPartitionPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {}
	virtual ~ParallelPartitionPartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		
		std::vector<std::set<int> > cons_double(model->population_size,std::set<int>());
		std::vector<idx_t> vwgt_v(model->population_size,1);
		
		int pre_size = model->population_size;
		for(int ii=0; ii < pre_size; ii++) {
			int from = ii;
			int post_size = model->interconnections_size[ii];
			for(int jj=0; jj < post_size; jj++) {
				int to = abs(model->interconnections[ii][jj]);
				// ignore looping connections (pre and post neurons are identical)
				// causes ParMETIS non-symmetric adjacency matrix error "key X not found"
				if(from == to) continue;
				// Take into consideration synaptic computational load to weight vertices load 
				vwgt_v[to]++;
				cons_double[from].insert(to);
				cons_double[to].insert(from);
			}
		}
		
		///////
		// create CSR representation of neurons and connections
		///////
		std::vector<idx_t> xadj_v(model->population_size + 1);
		std::vector<idx_t> adjncy_v;
		std::vector<idx_t> adjwgt_v;
		// create data structures for partitioning
		idx_t connection_index = 0;
		pre_size = cons_double.size();
		for(int ii=0; ii < pre_size; ii++) {
			xadj_v[ii] = connection_index; // WHAT HAPPENS IF A NEURON DOES NOT HAVE ANY OUTWARD CONNECTION??
			int post_size = cons_double[ii].size();
			if(post_size == 0) {
				adjncy_v.push_back(ii);
				adjwgt_v.push_back(1);
				connection_index++;
			} else {
				for (std::set<int>::iterator it = cons_double[ii].begin(); it != cons_double[ii].end(); ++it) {
					adjncy_v.push_back(*it);
					// check weight of connection ii -> cons_double[ii][jj]
					//int weight = 0; // default weight value for reciprocal connections (METIS requires bidirectional edges)
					//for (int dd=0; dd < model->connections[ii].size(); dd++) {
					//	if(model->connections[ii][dd] == cons_double[ii][jj]) {
					//		weight += abs(model->connection_weights[ii][dd]);
					//		break;
					//	}
					//}
					//for (int kk=0; kk < model->connections[cons_double[ii][jj]].size(); kk++) {
					//	if(model->connections[cons_double[ii][jj]][kk] == ii) {
					//		weight += abs(model->connection_weights[cons_double[ii][jj]][kk]);
					//		break;
					//	}
					//}
					//adjwgt_v.push_back(weight);
					// approximate the weight of each edge as the neuron activity from previous sims
					//adjwgt_v.push_back(previous_activity->at(ii)); // why is the sum needed?? otherwise segfault !?
					int weight = max(previous_activity->at(ii),1);
					adjwgt_v.push_back(weight); 
					connection_index++;
				}
			}
			
		}
		xadj_v[model->population_size] = connection_index;
		
		MPI_Comm comm = MPI_COMM_WORLD;
		idx_t nParts = partitions; // number of processes (partitions)
		idx_t ncon = 1;
		idx_t objval;
		
		//////
		// initially distribute vertices amongs processes
		//////
		int vert_per_process = floor((xadj_v.size() - 1) * 1.0f / partitions);
		std::vector<idx_t> vtxdist_v(partitions + 1);
		std::vector<idx_t> localxadj_v;
		std::vector<idx_t> localvwgt_v;
		int localVerts = 0;
		int minEdgeIndex = 0;
		int maxEdgeIndex = 0;
		for(int ii=0; ii < partitions; ii++) {
			vtxdist_v[ii] = localVerts;
			if(ii == process_id) {
				int baseIndex = xadj_v[localVerts];
				minEdgeIndex = baseIndex;
				//int maxIndex = localVerts + vert_per_process < xadj_v.size() ? vert_per_process : xadj_v.size() - localVerts - 1;
				int maxIndex = (ii == partitions-1) ? xadj_v.size() - localVerts - 1 : vert_per_process;
				for(int jj=0; jj < maxIndex; jj++) {
					localxadj_v.push_back(xadj_v[localVerts]-baseIndex);
					localvwgt_v.push_back(vwgt_v[localVerts]);
					localVerts++;
				}
				localxadj_v.push_back(xadj_v[localVerts]-baseIndex);
				maxEdgeIndex = xadj_v[localVerts];
			} else {
				localVerts += vert_per_process;
			}
		}
		vtxdist_v[partitions] = xadj_v.size() - 1;
		
		int total_adjncy = adjncy_v.size();
		if(maxEdgeIndex < adjncy_v.size()-1) {
			adjncy_v.erase(adjncy_v.end() - (adjncy_v.size()-maxEdgeIndex), adjncy_v.end());
			adjwgt_v.erase(adjwgt_v.end() - (adjwgt_v.size()-maxEdgeIndex), adjwgt_v.end());
		}
		if(minEdgeIndex > 0) {
			adjncy_v.erase(adjncy_v.begin(),adjncy_v.begin()+minEdgeIndex);
			adjwgt_v.erase(adjwgt_v.begin(),adjwgt_v.begin()+minEdgeIndex);
		}
		
		idx_t* xadj = &localxadj_v[0];
		idx_t* adjncy = &adjncy_v[0];
		idx_t* vwgt = &localvwgt_v[0];
		idx_t* vtxdist = &vtxdist_v[0];
		idx_t* adjwgt = &adjwgt_v[0];
		real_t* tpwgts;
		tpwgts = (real_t*)malloc(partitions * sizeof(real_t));
		for(int ii=0; ii < partitions; ii++) {
			tpwgts[ii] = 1.0 / partitions;
		}	
		
		if(nParts == 1) {
			PRINTF("Partitioning not required\n");
			for(int ii=0; ii < model->population_size; ii++) 
				partitioning[ii] = 0;
		} else {
			// partition using: vertices and edge weights
			// ParMETIS options
			idx_t numflag = 0; // 0 for C-style indices; 1 for fortran
			idx_t wgtflag = 3; // 3 if wgt in both vertices and edges
			real_t ubvec[1] = {1.003f}; // tolerance imbalance for each vertex weight
			idx_t options[3] = {
				1, 		// 0 for default options, 1 to read on
				0, 		// info returned from methods (see parmetis.h for values)
				rand()	// random seed				
			};
			
			PRINTF("%i: Calling ParMETIS\n",process_id);
			idx_t* local_partitioning = (idx_t*)calloc(localxadj_v.size()-1,sizeof(idx_t));
			int ret = ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, 
							&wgtflag, &numflag, &ncon, &nParts, 
							tpwgts, ubvec, options, &objval, local_partitioning, &comm);
		
			PRINTF("%i: objective value: %li\n",process_id,objval);
			
			// Now local_partitioning must be shared across processes
			int num_sends = localxadj_v.size()-1;
			int* recvcounts = (int*)malloc(sizeof(int) * partitions);
			int* displs = (int*)malloc(sizeof(int) * partitions);
			MPI_Allgather(&num_sends,1,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);
			int counter = 0;
			for(int ii=0; ii < partitions; ii++) {
				displs[ii] = counter;
				counter += recvcounts[ii]; 
			}
			MPI_Allgatherv(local_partitioning,num_sends,MPI_LONG,partitioning,recvcounts,displs,MPI_LONG,MPI_COMM_WORLD);
			free(local_partitioning);
		}
		
		// partitioning done, free resources
		free(tpwgts);			
	}
};


#endif
