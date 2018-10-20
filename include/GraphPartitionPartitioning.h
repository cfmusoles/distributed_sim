#ifndef GRAPH_PARTITION_PARTITIONING
#define GRAPH_PARTITION_PARTITIONING

#include <metis.h>
#include <vector>
#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"

class GraphPartitionPartitioning : public Partitioning {
public:
			
	GraphPartitionPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {}
	virtual ~GraphPartitionPartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		
		std::vector<idx_t> vsize_v(previous_activity->size());
		for(int ii=0; ii < previous_activity->size(); ii++) {
			vsize_v[ii] = (*previous_activity)[ii]; 
		}
		std::vector<std::vector<int> > cons_double(model->population_size,std::vector<int>(0));
		std::vector<idx_t> vwgt_v(model->population_size,1);
		
		int pre_size = model->population_size;
		for(int ii=0; ii < pre_size; ii++) {
			int from = ii;
			std::vector<int>* vec_from = &cons_double[from];
			int post_size = model->interconnections_size[ii];
			for(int jj=0; jj < post_size; jj++) {
				int to = abs(model->interconnections[ii][jj]);
				// Take into consideration synaptic computational load to weight vertices load 
				vwgt_v[to]++;
				// ignore looping connections (pre and post neurons are identical)
				// causes ParMETIS non-symmetric adjacency matrix error "key X not found"
				// METIS works fine, but for fairness we will not consider them (they should not affect partitioning)
				if(from == to) continue;
				std::vector<int>* vec_to = &cons_double[to];
				
				if(std::find(vec_from->begin(),vec_from->end(),to) == vec_from->end()) {
					vec_from->push_back(to);
				}
				if(std::find(vec_to->begin(),vec_to->end(),from) == vec_to->end()) {
					vec_to->push_back(from);
				}
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
				for(int jj=0; jj < post_size; jj++) {
					adjncy_v.push_back(cons_double[ii][jj]);
					// check weight of connection ii -> cons_double[ii][jj]
					/*int weight = 0; // default weight value for reciprocal connections (METIS requires bidirectional edges)
					for (int dd=0; dd < model->connections[ii].size(); dd++) {
						if(model->connections[ii][dd] == cons_double[ii][jj]) {
							weight += abs(model->connection_weights[ii][dd]);
							break;
						}
					}
					for (int kk=0; kk < model->connections[cons_double[ii][jj]].size(); kk++) {
						if(model->connections[cons_double[ii][jj]][kk] == ii) {
							weight += abs(model->connection_weights[cons_double[ii][jj]][kk]);
							break;
						}
					}
					adjwgt_v.push_back(weight);*/
					adjwgt_v.push_back(1); // don't consider edge weights at the moment
					connection_index++;
				}
			}
			
		}
		xadj_v[model->population_size] = connection_index;
		
		//////
		// partition xadj and divide workload amongst MPI processes
		//////
		idx_t nParts = partitions; // number of processes (partitions)
		idx_t nWeights = 1;
		idx_t ncon = 1;
		idx_t objval;
		idx_t* xadj = &xadj_v[0];
		idx_t* adjncy = &adjncy_v[0];
		idx_t* vwgt = &vwgt_v[0];
		idx_t* vsize = &vsize_v[0];
		idx_t* adjwgt = &adjwgt_v[0];
		real_t* tpwgts;
		tpwgts = (real_t*)malloc(nParts * sizeof(real_t));
		for(int ii=0; ii < nParts; ii++) {
			tpwgts[ii] = 1.0 / nParts;
		}	
		
		
		if(nParts == 1) {
			PRINTF("Partitioning not required\n");
			for(int ii=0; ii < model->population_size; ii++) 
				partitioning[ii] = 0;
		} else {
			// partition using: vertices and edge weights
			// METIS options
			idx_t options[METIS_NOPTIONS];
			options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // METIS_PTYPE_KWAY, METIS_PTYPE_RB 
			options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // METIS_OBJTYPE_VOL, METIS_OBJTYPE_CUT
			options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // METIS_CTYPE_SHEM, METIS_CTYPE_RM
			options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE; // METIS_IPTYPE_EDGE, METIS_IPTYPE_GROW, METIS_IPTYPE_RANDOM, METIS_IPTYPE_NODE
			options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // METIS_RTYPE_FM, METIS_RTYPE_GREEDY, METIS_RTYPE_SEP2SIDED, METIS_RTYPE_SEP1SIDED
			options[METIS_OPTION_NCUTS] = 1;	// number of different partitionings performed (chosen one is the one with lowest edgecut)
			options[METIS_OPTION_NSEPS] = 1;	// number of separators at each level of nested dissection
			options[METIS_OPTION_NUMBERING] = 0; // 0 - Cstyle (starts from 0), 1 - Fortran (starts from 1)
			options[METIS_OPTION_NITER] = 100;		// number of iterations for refinement algorithm during uncoarsening
			options[METIS_OPTION_SEED] = rand();
			options[METIS_OPTION_MINCONN] = 0; // 0 - not minimise maximum connectivity, 1- minimise maximum connectivity
			options[METIS_OPTION_NO2HOP] = 0; // 0- Use 2-hop matching, 1 - not use it
			options[METIS_OPTION_CONTIG] = 0; // 0 - Don't force contiguous partitions, 1 - force contiguity
			options[METIS_OPTION_COMPRESS] = 0; // 0 - No compression, 1 - compression
			options[METIS_OPTION_PFACTOR] = 0;
			options[METIS_OPTION_UFACTOR] = 3; // load imbalance factor (1 + x) / 1000
			options[METIS_OPTION_DBGLVL] = 0;
			
				
			PRINTF("%i: Calling METIS\n",process_id);
			int ret = METIS_PartGraphKway(&(model->population_size), &ncon, xadj, adjncy,
							vwgt, vsize, adjwgt, &nParts, tpwgts,
							NULL, options, &objval, partitioning);
		
			PRINTF("%i: objective value: %li\n",process_id,objval);
		}
		
		// partitioning done, free resources
		free(tpwgts);			
	}
};


/*
class GraphPartitionPartitioning : public Partitioning {
public:
	std::vector<std::vector<int> > cons_double;
	std::vector<idx_t> vsize_v;
	std::vector<idx_t> vwgt_v;
			
	GraphPartitionPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {}
	virtual ~GraphPartitionPartitioning() {}
	
	virtual void clear_intermediate_data()  {
		std::vector<std::vector<int> >().swap(cons_double);
		std::vector<idx_t>().swap(vsize_v);
		std::vector<idx_t>().swap(vwgt_v);
	}
	
	virtual void generate_connection_structures(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		vsize_v.clear();
		vsize_v.resize(previous_activity->size());
		for(int ii=0; ii < previous_activity->size(); ii++) {
			vsize_v[ii] = ((*previous_activity)[ii]); 
		}
		cons_double.resize(model->population_size);
		vwgt_v.clear();
		vwgt_v.resize(model->population_size,1);
		
		for(int ii=0; ii < model->connections.size(); ii++) {
			int from = ii;
			for(int jj=0; jj < model->connections[ii].size(); jj++) {
				int to = model->connections[ii][jj];
				// Take into consideration synaptic computational load to weight vertices load 
				vwgt_v[to]++;
				if(std::find(cons_double[from].begin(),cons_double[from].end(),to) == cons_double[from].end()) {
					cons_double[from].push_back(to);
				}
				if(std::find(cons_double[to].begin(),cons_double[to].end(),from) == cons_double[to].end()) {
					cons_double[to].push_back(from);
				}
			}
		}
		
	}
	
	virtual void assign_partition_units(Model* model, int partitions, int process_id) {
		
		///////
		// create CSR representation of neurons and connections
		///////
		std::vector<idx_t> xadj_v(model->population_size + 1);
		std::vector<idx_t> adjncy_v;
		std::vector<idx_t> adjwgt_v;
		// create data structures for partitioning
		idx_t connection_index = 0;
		for(int ii=0; ii < cons_double.size(); ii++) {
			xadj_v[ii] = connection_index; // WHAT HAPPENS IF A NEURON DOES NOT HAVE ANY OUTWARD CONNECTION??
			if(cons_double[ii].size() == 0) {
				adjncy_v.push_back(ii);
				adjwgt_v.push_back(1);
				connection_index++;
			} else {
				for(int jj=0; jj < cons_double[ii].size(); jj++) {
					adjncy_v.push_back(cons_double[ii][jj]);
					// check weight of connection ii -> cons_double[ii][jj]
					int weight = 0; // default weight value for reciprocal connections (METIS requires bidirectional edges)
					for (int dd=0; dd < model->connections[ii].size(); dd++) {
						if(model->connections[ii][dd] == cons_double[ii][jj]) {
							weight += abs(model->connection_weights[ii][dd]);
							break;
						}
					}
					for (int kk=0; kk < model->connections[cons_double[ii][jj]].size(); kk++) {
						if(model->connections[cons_double[ii][jj]][kk] == ii) {
							weight += abs(model->connection_weights[cons_double[ii][jj]][kk]);
							break;
						}
					}
					//adjwgt_v.push_back(weight);
					adjwgt_v.push_back(1); // don't consider edge weights at the moment
					connection_index++;
				}
			}
			
		}
		xadj_v[model->population_size] = connection_index;
		
		//////
		// partition xadj and divide workload amongst MPI processes
		//////
		idx_t nParts = partitions; // number of processes (partitions)
		idx_t nWeights = 1;
		idx_t ncon = 1;
		idx_t objval;
		idx_t* xadj = &xadj_v[0];
		idx_t* adjncy = &adjncy_v[0];
		idx_t* vwgt = &vwgt_v[0];
		idx_t* vsize = &vsize_v[0];
		idx_t* adjwgt = &adjwgt_v[0];
		real_t* tpwgts;
		tpwgts = (real_t*)malloc(nParts * sizeof(real_t));
		for(int ii=0; ii < nParts; ii++) {
			tpwgts[ii] = 1.0 / nParts;
		}	
		
		
		if(nParts == 1) {
			PRINTF("Partitioning not required\n");
			for(int ii=0; ii < model->population_size; ii++) 
				partitioning[ii] = 0;
		} else {
			// partition using: vertices and edge weights
			// METIS options
			idx_t options[METIS_NOPTIONS];
			options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // METIS_PTYPE_KWAY, METIS_PTYPE_RB 
			options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // METIS_OBJTYPE_VOL, METIS_OBJTYPE_CUT
			options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // METIS_CTYPE_SHEM, METIS_CTYPE_RM
			options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE; // METIS_IPTYPE_EDGE, METIS_IPTYPE_GROW, METIS_IPTYPE_RANDOM, METIS_IPTYPE_NODE
			options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // METIS_RTYPE_FM, METIS_RTYPE_GREEDY, METIS_RTYPE_SEP2SIDED, METIS_RTYPE_SEP1SIDED
			options[METIS_OPTION_NCUTS] = 10;	// number of different partitionings performed (chosen one is the one with lowest edgecut)
			options[METIS_OPTION_NSEPS] = 1;	// number of separators at each level of nested dissection
			options[METIS_OPTION_NUMBERING] = 0; // 0 - Cstyle (starts from 0), 1 - Fortran (starts from 1)
			options[METIS_OPTION_NITER] = 100;		// number of iterations for refinement algorithm during uncoarsening
			options[METIS_OPTION_SEED] = rand();
			options[METIS_OPTION_MINCONN] = 0; // 0 - not minimise maximum connectivity, 1- minimise maximum connectivity
			options[METIS_OPTION_NO2HOP] = 0; // 0- Use 2-hop matching, 1 - not use it
			options[METIS_OPTION_CONTIG] = 0; // 0 - Don't force contiguous partitions, 1 - force contiguity
			options[METIS_OPTION_COMPRESS] = 0; // 0 - No compression, 1 - compression
			options[METIS_OPTION_PFACTOR] = 0;
			options[METIS_OPTION_UFACTOR] = 5; // load imbalance factor (1 + x) / 1000
			options[METIS_OPTION_DBGLVL] = 0;
			
				
			PRINTF("%i: Calling METIS\n",process_id);
			//int ret = METIS_PartGraphKway(&(model->population_size), &ncon, xadj, adjncy,
			//				vwgt, vsize, adjwgt, &nParts, tpwgts,
			//				NULL, options, &objval, partitioning);
		
			PRINTF("%i: objective value: %li\n",process_id,objval);
		}
		
		// partitioning done, free resources
		free(tpwgts);
		
	}
};*/

#endif
