#ifndef ZOLTAN_HYPERGRAPH_PARTITION_PARTITIONING
#define ZOLTAN_HYPERGRAPH_PARTITION_PARTITIONING

#include "zoltan.h"
#include <vector>
#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"
#include "HypergraphPartitioning.h"
#include "PRAW.h"
#include <cstring>
#include <cstdio>
#include <stdlib.h>

class ZoltanFilePartitioning : public Partitioning {
public:
	
	ZoltanFilePartitioning(std::vector<Population*>* pops, int population_size,char* comm_bandwidth_file) : Partitioning(pops,population_size) {
		comm_bandwidth_filename = comm_bandwidth_file;
        float ver;
		rc = Zoltan_Initialize(0, NULL, &ver);
	
		if (rc != ZOLTAN_OK){
			printf("Error initialising Zoltan...\n");
			MPI_Finalize();
			exit(0);
		}

		zz = Zoltan_Create(MPI_COMM_WORLD);

		setZoltanParams();
		// initialise hg struct
		hg.vtxGID = NULL;
		hg.vtx_wgts = NULL;
		hg.vtxedge_ptr = NULL;
		hg.pin_GID = NULL;
	}
	virtual ~ZoltanFilePartitioning() {
		Zoltan_Destroy(&zz);
	}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		
		// Assign unique vertices to partitions
		hg.numMyVertices = 0;
		for(int ii=0; ii < model->population_size; ii++) {
            partitioning[ii] = (double)rand() / (double)RAND_MAX * partitions;
            if(partitioning[ii] == partitions) partitioning[ii] -= 1;
			if(partitioning[ii] == process_id) hg.numMyVertices++;
        }
        
        if(model->hypergraph_file == NULL) {
            PRINTF("%i: hypergraph file not set in Model object. Random partition.",process_id);
            return;
        }

		// 1 Load model from hmetis file model->hypergraph_file
		std::vector<std::vector<int> > hyperedges(model->population_size);
        std::vector<std::vector<int> > hedge_ptr(model->population_size);
        PRAW::load_hypergraph_from_file(model->hypergraph_file,&hyperedges,&hedge_ptr);

		// 2 Convert hmetis format to zoltan
		hg.process_id = process_id; 
		hg.partitioning = partitioning;

		int* total_vtx_wgts_v = (int*)calloc(model->population_size,sizeof(int));
		hg.vtxGID = (ZOLTAN_ID_TYPE*)malloc(sizeof(ZOLTAN_ID_TYPE) * hg.numMyVertices);
		hg.vtxedge_ptr = (int*)malloc(sizeof(int) * (hg.numMyVertices+1));
		int current_local_vertex = 0;
		hg.npins = 0;
		for(int from=0; from < model->population_size; from++) {
			total_vtx_wgts_v[from] += 1; // accounts for neuron update computation 
			if(partitioning[from] == process_id) {
				hg.vtxGID[current_local_vertex] = from;
				hg.vtxedge_ptr[current_local_vertex] = hg.npins;
				current_local_vertex++;
				for(int jj=0; jj < hedge_ptr[from].size(); jj++) {
					int he_id = hedge_ptr[from][jj];
					for(int pp=0; pp < hyperedges[he_id].size(); pp++) {
						int to = hyperedges[he_id][pp];
                        if(to == from) continue; //
						hg.npins++;
						// Take into consideration synaptic computational load to weight vertices load 
						// issue with hmetis hypergraph --> a connection between two neurons may exist more than once in the model->intereconnections array
						// incoming synaptic computational balance must only count each once
						total_vtx_wgts_v[to] += hyperedges[he_id].size()-1;
					}
					
				}
				
			}
		}
		hg.vtxedge_ptr[hg.numMyVertices] = hg.npins;
		hg.pin_GID = (ZOLTAN_ID_TYPE*)malloc(sizeof(ZOLTAN_ID_TYPE) * hg.npins);
		int counter = 0;
		for(int from=0; from < model->population_size; from++) {
			if(partitioning[from] == process_id) {
				for(int jj=0; jj < hedge_ptr[from].size(); jj++) {
					int he_id = hedge_ptr[from][jj];
					for(int pp=0; pp < hyperedges[he_id].size(); pp++) {
						if(from == hyperedges[he_id][pp]) continue; //
						hg.pin_GID[counter] = hyperedges[he_id][pp];
						counter++;
					}
				}
				
			}
		}
		
		hg.vtx_wgts = (float*)malloc(sizeof(float) * hg.numMyVertices);
		current_local_vertex = 0;
		for(int ii=0; ii < model->population_size; ii++) {
			if(partitioning[ii] == process_id) {
				//hg.vtx_wgts[current_local_vertex] = 1;
				hg.vtx_wgts[current_local_vertex] = model->null_compute ? 1 : total_vtx_wgts_v[ii];
				current_local_vertex++;
			}
		}
		free(total_vtx_wgts_v);
		
		

		// 3 Do partitioning
		int changes, numGidEntries, numLidEntries, numImport, numExport;
		ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
		int *importProcs, *importToPart, *exportProcs, *exportToPart;
		
		PRINTF("%i: Calling Zoltan\n",process_id);
		// Zoltan_LB_Set_Part_Sizes to specify the size of the partitions http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_lb.html#Zoltan_LB_Set_Part_Sizes
		// Zoltan_LB_Eval to evaluate quality of decomposition http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_lb.html#Zoltan_LB_Eval
		//ZOLTAN_BALANCE_EVAL balance_eval;
		//ZOLTAN_GRAPH_EVAL graph_eval;
		//ZOLTAN_HG_EVAL hg_eval;
		//Zoltan_LB_Eval(zz,1,&balance_eval,&graph_eval,&hg_eval);
			
		rc = Zoltan_LB_Partition(zz, // input (all remaining fields are output) 
			&changes,        // 1 if partitioning was changed, 0 otherwise  
			&numGidEntries,  // Number of integers used for a global ID 
			&numLidEntries,  // Number of integers used for a local ID 
			&numImport,      // Number of vertices to be sent to me 
			&importGlobalGids,  // Global IDs of vertices to be sent to me 
			&importLocalGids,   // Local IDs of vertices to be sent to me 
			&importProcs,    // Process rank for source of each incoming vertex 
			&importToPart,   // New partition for each incoming vertex 
			&numExport,      // Number of vertices I must send to other processes
			&exportGlobalGids,  // Global IDs of the vertices I must send 
			&exportLocalGids,   // Local IDs of the vertices I must send 
			&exportProcs,    // Process to which I send each of the vertices 
			&exportToPart);  // Partition to which each vertex will belong 
		
		if (rc != ZOLTAN_OK){
			printf("Zoltan error after partitioning...\n");
			MPI_Finalize();
			Zoltan_Destroy(&zz);
			exit(0);
		}
		
		// We are doing AUTO_MIGRATE
		// http://www.cs.sandia.gov/zoltan/ug_html/ug_interface_mig.html
		// after migration, partitioning within each part needs to be purged for non-owned vertices 
		// as the information may no longer be accurate
		
		idx_t* partAssign = (idx_t*)calloc(model->population_size,sizeof(idx_t)); // calloc will init to 0 all
		for (int i=0; i < model->population_size; i++) {
			if(partitioning[i] == process_id)
				partAssign[i] = process_id;		
		}
		MPI_Allreduce(partAssign,partitioning,model->population_size,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
		free(partAssign);

		if(process_id == 0) {
            std::string filename = model->hypergraph_file;
			filename += "_zoltanFile";
			// p2p communication cost estimates from file
			double** comm_cost_matrix = (double**)malloc(sizeof(double*) * partitions);
			for(int ii=0; ii < partitions; ii++) {
				comm_cost_matrix[ii] = (double*)calloc(partitions,sizeof(double));
			}
			if(comm_bandwidth_filename != NULL) {
				PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,partitions,false);
			}
			PRAW::storePartitionStats(filename,partitioning,partitions,model->population_size,&hyperedges,&hedge_ptr,NULL,comm_cost_matrix);
			
			for(int ii=0; ii < partitions; ii++) {
				free(comm_cost_matrix[ii]);
			}
			free(comm_cost_matrix);
		}
		Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
						&importProcs, &importToPart);
		Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
						&exportProcs, &exportToPart);	
		if(hg.vtxGID != NULL) free(hg.vtxGID);
		if(hg.pin_GID != NULL) free(hg.pin_GID);
		if(hg.vtxedge_ptr != NULL) free(hg.vtxedge_ptr);
		if(hg.vtx_wgts != NULL) free(hg.vtx_wgts);

		
	}


private:

	struct Zoltan_Struct *zz;	
	HGRAPH_DATA hg;
	int rc;
	char* comm_bandwidth_filename = NULL;

	void setZoltanParams() {
		/* General parameters */
	
		Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
		Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
		Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
		Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
		Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
		Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
		Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* use Zoltan default vertex weights */
		Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */
		Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.1");/* imbalance tolerance */

		// parameters for PHG: http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_phg.html
		
		Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
		Zoltan_Set_Param(zz, "CHECK_HYPERGRAPH", "0");
		Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "MAX");
		Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", "1.0");
		Zoltan_Set_Param(zz, "PHG_CUT_OBJECTIVE", "HYPEREDGES"); // CONNECTIVITY, HYPEREDGES
		Zoltan_Set_Param(zz, "PHG_COARSENING_METHOD", "AGG");
		Zoltan_Set_Param(zz, "PHG_COARSEPARTITION_METHOD", "AUTO");
		Zoltan_Set_Param(zz, "PHG_REFINEMENT_METHOD", "FM");
		Zoltan_Set_Param(zz, "PHG_REFINEMENT_QUALITY", "1");
		Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");  
		
		/* Application defined query functions */
		/* To set partitioning callbacks, see: http://www.cs.sandia.gov/zoltan/ug_html/ug_query_lb.html */
		Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &hg);	// number of objects owned by processor
		Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &hg);			// list of weights for owned objects
		Zoltan_Set_HG_Size_CS_Fn(zz, get_hypergraph_size, &hg);	// dimensionality of the problem
		Zoltan_Set_HG_CS_Fn(zz, get_hypergraph, &hg);				//  coordinate system (or connectivity)
		//Zoltan_Set_HG_Edge_Wts_Fn(zz, get_hyperedges_weights, &hg);  // get hyperedges weights
		//Zoltan_Set_HG_Size_Edge_Wts_Fn(zz, get_hyperedges_size_weights, &hg); // get length of hyperedges weights
		/* for migration */
		Zoltan_Set_Obj_Size_Multi_Fn(zz, user_migration_multi_object_size, &hg);
		Zoltan_Set_Pack_Obj_Multi_Fn(zz, user_pack_multi_obj, &hg);
		Zoltan_Set_Unpack_Obj_Multi_Fn(zz, user_unpack_multi_obj, &hg);
	}
};


#endif
