/* Hierarchical Hypergraph partitioning class
* Info on: http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_hier.html
*/

#ifndef HIER_PARTITION_PARTITIONING
#define HIER_PARTITION_PARTITIONING

#include "zoltan.h"
#include <vector>
#include <set>
#include "Simulation.h"
#include "HypergraphPartitioning.h"
#include "Utils.h"

//////////////////////////////////////////////////
// ZOLTAN SPECIFIC DATA STRUCTURES AND FUNCTIONS for Hierarchical partitioning
// General Hypergraph functions are defined in HypergraphPartitioning
//////////////////////////////////////////////////

/* Application defined query functions. */
#pragma region Zoltan_app_functions

// for the calling processor, the number of levels of hierarchy for hierarchical load balancing
int hier_num_levels_fn(void* data, int* ierr) {
    // called only once
	*ierr = ZOLTAN_OK;
    return 2;
}

// gets the part number to be used for the given level of a hierarchical balancing procedure
int hier_part_fn(void* data, int level, int* ierr) {
    // at each level, partition number to be computed at that level for calling process
    // called once per level
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	if(hg == NULL) {
		*ierr = ZOLTAN_FATAL;
		printf("hier_part_fn: HGRAPH_DATA passed was null");
		return -1;
	} 
	
	switch(level) {
		case 0:
		{
			return hg->process_id / 2;
		}
		case 1:
		{
			return hg->process_id % 2;
		}
	}
    *ierr = ZOLTAN_OK;
}

//provides to the calling process the Zoltan_Struct to be used to guide the partitioning
// and load balancing at the given level in the hierarchy. 
// This Zoltan_Struct can be passed to Zoltan_Set_Param to set load balancing parameters 
// for this level in the hierarchical balancing
void hier_method_fn(void* data, int level, Zoltan_Struct* zz, int* ierr) {
    // called once per level. Can be different method at each level
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
	Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
	Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* use Zoltan default vertex weights */
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.001");/* imbalance tolerance */

	// parameters for PHG: http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_phg.html
	Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
	Zoltan_Set_Param(zz, "CHECK_HYPERGRAPH", "1");
	Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "MAX");
	Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", "1.0");
	Zoltan_Set_Param(zz, "PHG_CUT_OBJECTIVE", "CONNECTIVITY"); // CONNECTIVITY, HYPEREDGES
	Zoltan_Set_Param(zz, "PHG_COARSENING_METHOD", "AGG");
	Zoltan_Set_Param(zz, "PHG_COARSEPARTITION_METHOD", "AUTO");
	Zoltan_Set_Param(zz, "PHG_REFINEMENT_METHOD", "FM");
	Zoltan_Set_Param(zz, "PHG_REFINEMENT_QUALITY", "10");

	// hg is correct in first level but not on second!
	Zoltan_Set_HG_Size_CS_Fn(zz, get_hypergraph_size, hg);	// dimensionality of the problem
	Zoltan_Set_HG_CS_Fn(zz, get_hypergraph, hg);				//  coordinate system (or connectivity)
	*ierr = ZOLTAN_OK;
}

#pragma endregion

class HIERPartitioning : public Partitioning {
public:
	HIERPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {
		float ver;
		rc = Zoltan_Initialize(0, NULL, &ver);
	
		if (rc != ZOLTAN_OK){
			printf("Error initialising Zoltan...\n");
			MPI_Finalize();
			exit(0);
		}

		zz = Zoltan_Create(MPI_COMM_WORLD);

		setZoltanParams();		
	}
	virtual ~HIERPartitioning() {
		Zoltan_Destroy(&zz);
	}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		
		hg.process_id = process_id; 
		hg.partitioning = partitioning;
		
		// Initially distributed assignment
		int current_assignment = 0;
		hg.numMyVertices = 0;
		
		for(int ii=0; ii < model->population_size; ii++) {
			if(current_assignment >= partitions) current_assignment = 0;
			partitioning[ii] = current_assignment;
			if(current_assignment == process_id) hg.numMyVertices++;
			current_assignment++;
		}
		
		if(partitions <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

		// generate local vtxedge_GID and global vtx weights and hyperedge ids
		int* total_vtx_wgts_v = (int*)calloc(model->population_size,sizeof(int));
		hg.vtxGID = (ZOLTAN_ID_TYPE*)malloc(sizeof(ZOLTAN_ID_TYPE) * hg.numMyVertices);
		int* hyperedges_sizes = (int*)calloc(model->population_size,sizeof(int));
		int pre_size = model->population_size;
		int post_size, to;
		int current_local_vertex = 0;
		hg.npins = 0;
		for(int from=0; from < pre_size; from++) {
			total_vtx_wgts_v[from] += 1; // accounts for neuron update computation 
			if(partitioning[from] == process_id) {
				hg.vtxGID[current_local_vertex] = from;
				current_local_vertex++;
				hg.npins++;
				hyperedges_sizes[from] += 1;
			}
			post_size = model->interconnections_size[from];
			for(int jj=0; jj < post_size; jj++) {
				to = abs(model->interconnections[from][jj]);
				// Take into consideration synaptic computational load to weight vertices load 
				total_vtx_wgts_v[to] += 1; // accounts for synaptic update computation
				if(from == to) continue;
				if(partitioning[to] == process_id) {
					hg.npins++; 
					hyperedges_sizes[to] += 1;
				}
			}
		}
		// create hyperedges
		int** hyperedges = (int**)malloc(sizeof(int*) * model->population_size);
		for(int ii=0; ii < model->population_size; ii++) {
			if(hyperedges_sizes[ii] > 0) hyperedges[ii] = (int*)malloc(sizeof(int) * hyperedges_sizes[ii]);
		}
		int* lastId = (int*)calloc(model->population_size,sizeof(int));
		for(int from=0; from < pre_size; from++) {
			if(partitioning[from] == process_id) {
				hyperedges[from][lastId[from]] = from;
				lastId[from] += 1;
			}
			post_size = model->interconnections_size[from];
			for(int jj=0; jj < post_size; jj++) {
				to = abs(model->interconnections[from][jj]);
				if(from == to) continue;
				if(partitioning[to] == process_id) {
					hyperedges[to][lastId[to]] = from;
					lastId[to] += 1;
				}
			}
			
		}
		free(lastId);
		
		// create local pin_GID and vtxedge_ptr vectors
		hg.pin_GID = (ZOLTAN_ID_TYPE*)malloc(sizeof(ZOLTAN_ID_TYPE) * hg.npins);
		hg.vtxedge_ptr = (int*)malloc(sizeof(int) * (hg.numMyVertices+1));
		int index = 0;
		current_local_vertex = 0;
		for(int from=0; from < pre_size; from++) {
			if(partitioning[from] == process_id) {
				hg.vtxedge_ptr[current_local_vertex] = index;
				post_size = hyperedges_sizes[from];
				for(int jj=0; jj < post_size; jj++) {
					hg.pin_GID[index] = hyperedges[from][jj];
					index++;
				}
				current_local_vertex++;
			}
			
		}
		hg.vtxedge_ptr[current_local_vertex] = index;
		
		// clean up intermediate data structures
		for(int ii=0; ii < model->population_size; ii++) {
			if(hyperedges_sizes[ii] > 0) free(hyperedges[ii]);
		}
		free(hyperedges);
		free(hyperedges_sizes);
		
		// refine total_vtx_wgts_v to only include local references
		hg.vtx_wgts = (float*)malloc(sizeof(float) * hg.numMyVertices);
		current_local_vertex = 0;
		for(int ii=0; ii < model->population_size; ii++) {
			if(partitioning[ii] == process_id) {
				hg.vtx_wgts[current_local_vertex] = total_vtx_wgts_v[ii];
				current_local_vertex++;
			}
		}
		free(total_vtx_wgts_v);
		
		//hg.numMyHEdges = vtxedge_GID_v.size(); // one hyperedge per presynaptic neuron
		//hg.edge_wgts = (float*)malloc(vtxedge_GID_v.size() * sizeof(float));
		//for(int ii=0; ii < vtxedge_GID_v.size(); ii++) 
		//	hg.edge_wgts = 1.0f;
		//hg.edgeGID // should be a list of unique hyperedges (from 0 to numMyHEdges; only local!)
		
		// partition using: Zoltan hypergraph
		int changes, numGidEntries, numLidEntries, numImport, numExport;
		ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
		int *importProcs, *importToPart, *exportProcs, *exportToPart;
		
		printf("%i: Prior partitioning value %lu\n",hg.process_id,hg.partitioning[19]);

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

		printf("%i: After partitioning value %lu\n",hg.process_id,hg.partitioning[19]);
		
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

	void setZoltanParams() {
		/// TODO: there seems to be a problem with registering callback functions here and
		// in different hier levels

		// General parameters 
        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
        Zoltan_Set_Param(zz, "HIER_DEBUG_LEVEL", "2");		
		Zoltan_Set_Param(zz, "LB_METHOD", "HIER");  	// partitioning method 
        Zoltan_Set_Param(zz, "HIER_ASSIST", "0");	// 1 for homogeneous multicore partitioning 
  		//Zoltan_Set_Param(zz, "TOPOLOGY", "2");  //"2,4,6" may refer to a dual-socket, 4-die, 6-core node 
		Zoltan_Set_Param(zz, "HIER_CHECKS", "1");   	// perform sanity checks	
		
		Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");  
		Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");// global IDs are integers 
        Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");// local IDs are integers 
        Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); // export AND import lists 
        Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); // use Zoltan default vertex weights 
        Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");// use Zoltan default hyperedge weights 
        Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.001");// imbalance tolerance 

		/* Application defined query functions */
		/* To set partitioning callbacks, see: http://www.cs.sandia.gov/zoltan/ug_html/ug_query_lb.html */
		Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &hg);	// number of objects owned by processor
		Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &hg);			// list of weights for owned objects
		// for migration 
		Zoltan_Set_Obj_Size_Multi_Fn(zz, user_migration_multi_object_size, &hg);
		Zoltan_Set_Pack_Obj_Multi_Fn(zz, user_pack_multi_obj, &hg);
		Zoltan_Set_Unpack_Obj_Multi_Fn(zz, user_unpack_multi_obj, &hg);
		// HIER partitioning methods
        Zoltan_Set_Hier_Num_Levels_Fn(zz, hier_num_levels_fn, &hg);
        Zoltan_Set_Hier_Part_Fn(zz, hier_part_fn, &hg);
        Zoltan_Set_Hier_Method_Fn(zz, hier_method_fn, &hg);

	}
};


#endif
