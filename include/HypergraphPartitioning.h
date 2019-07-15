#ifndef HYPERGRAPH_PARTITION_PARTITIONING
#define HYPERGRAPH_PARTITION_PARTITIONING

#include "zoltan.h"
#include <vector>
#include <algorithm>
#include <set>
#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"
#include <cstring>
#include <cstdio>
#include <stdlib.h>

//////////////////////////////////////////////////
// ZOLTAN SPECIFIC DATA STRUCTURES AND FUNCTIONS
//////////////////////////////////////////////////

typedef struct{
	/* Zoltan will partition vertices, while minimizing edge cuts */
	int numMyVertices;  /* number of vertices that I own initially */
	ZOLTAN_ID_TYPE *vtxGID;        /* global ID of these vertices */
	float* vtx_wgts;	/* weights of each of the vertices */
	//int numMyHEdges;    /* number of my hyperedges */
	//float* edge_wgts;	/* weight of hyperedges */
	int npins; /* number of vertices in my hyperedges */
	//ZOLTAN_ID_TYPE *edgeGID;       /* global ID of each of my hyperedges */
	int *vtxedge_ptr;     /* index into pin_GID array of edge's vertices */ // vtxedge_ptr
	ZOLTAN_ID_TYPE *pin_GID;  /* Edges of vertex vtxGID[i] begin at pin_GID[vtxedge_ptr[i]] */ // pin_GID
	int process_id;	// what process this graph resides in
	idx_t* partitioning;	// pointer to the partitioning array
} HGRAPH_DATA;

/* Application defined query functions. */
#pragma region Zoltan_app_functions
int get_number_of_vertices(void *data, int *ierr) {
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	PRINTF("%i: get_number_of_vertices\n\n",hg->process_id);
	*ierr = ZOLTAN_OK;
	return hg->numMyVertices;
}
void get_vertex_list(void *data, int num_gid_entries, int num_lid_entries,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr) {
	HGRAPH_DATA *hg= (HGRAPH_DATA *)data;
	PRINTF("%i: get_vertex_list\n\n",hg->process_id);

	*ierr = ZOLTAN_OK;
	/* For weights: 
	* 	wgt_dim is the number of weights (set by parameter OBJ_WEIGHT_DIM)
	* 	obj_wgts: associated weights. Weights for object i are stored in obj_wgts[(i-1)*wgt_dim:i*wgt_dim-1]
	*/
	for (int i=0; i<hg->numMyVertices; i++){
		globalID[i] = hg->vtxGID[i];
		localID[i] = i;
		if(wgt_dim > 0) obj_wgts[i] = hg->vtx_wgts[i];
	}

}
void get_hypergraph_size(void *data, int *num_lists, int *num_pins,
                                int *format, int *ierr){
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	PRINTF("%i: get_hypergraph_size (%i, %i)\n\n",hg->process_id,hg->numMyVertices,hg->npins);
	*ierr = ZOLTAN_OK;
	*num_lists = hg->numMyVertices; // number of vertices
	*num_pins = hg->npins; // number of elements in pin_GID
	*format = ZOLTAN_COMPRESSED_VERTEX;
	
}
void get_hypergraph(void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
                           int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
                           ZOLTAN_ID_PTR pin_GIDs, int *ierr) {
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	PRINTF("%i: get_hypergraph\n\n",hg->process_id);
	*ierr = ZOLTAN_OK;

	if ( (num_vtx_edge != hg->numMyVertices) || (num_pins != hg->npins) ||
	   (format != ZOLTAN_COMPRESSED_VERTEX)) {
		*ierr = ZOLTAN_FATAL;
		return;
	}
	int i;
	for (i=0; i < num_vtx_edge; i++){
		vtxedge_GID[i] = hg->vtxGID[i];
		vtxedge_ptr[i] = hg->vtxedge_ptr[i];
	}
	if(num_vtx_edge > 0)
		vtxedge_ptr[num_vtx_edge] = hg->vtxedge_ptr[num_vtx_edge];
	for (i=0; i < num_pins; i++){
		pin_GIDs[i] = hg->pin_GID[i];
	}
}
/*void get_hyperedges_weights(void *data, int num_gid_entries, int num_lid_entries, 
							int num_edges, int edge_weight_dim, ZOLTAN_ID_PTR edge_GID, 
							ZOLTAN_ID_PTR edge_LID, float  *edge_weight, int *ierr) {
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	*ierr = ZOLTAN_OK;
	// if each processor only sends weights for subset of edges, it does not need to be equal
	// actually, since numMyHEdges contains only the processor edges, this may be alright!
	if ( (num_edges != hg->numMyHEdges)) { 
		*ierr = ZOLTAN_FATAL;
		return;
	}
	// need to set edge_GID, edge_LID and edge_weight 
	// http://www.cs.sandia.gov/zoltan/ug_html/ug_query_lb.html#ZOLTAN_HG_EDGE_WTS_FN
	for (int i=0; i < num_edges; i++){ 
		edge_GID[i] = hg->edgeGID[i];
		edge_LID[i] = i;
		if(edge_weight_dim > 0) edge_weight[i] = hg->edge_wgts[i];
	} 
		
}
void get_hyperedges_size_weights(void *data, int *num_edges, int *ierr) {
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	*ierr = ZOLTAN_OK;
	*num_edges = hg->numMyHEdges;
}*/
void user_migration_multi_object_size(void *data, int num_gid_entries, int num_lid_entries, int num_ids, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,int* sizes, int *ierr) {
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	PRINTF("%i: user_migration_multi_object_size\n\n",hg->process_id);
	// set sizes. We will transfer the pins (connections) for each vertex and the vertex weight (represented as int here)
	for(int ii=0; ii < num_ids; ii++) {
		ZOLTAN_ID_TYPE local_id = local_ids[ii];
		int pins_for_vertex = hg->vtxedge_ptr[local_id+1] - hg->vtxedge_ptr[local_id];
		sizes[ii] = pins_for_vertex * sizeof(ZOLTAN_ID_TYPE) + sizeof(int);
	}
	*ierr = ZOLTAN_OK;
	
}
void user_pack_multi_obj(void *data, int num_gid_entries, int num_lid_entries, int num_ids, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int* dest, int* sizes, int* idx, char *buf, int *ierr) {
	
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	PRINTF("%i: user_pack_multi_obj\n\n",hg->process_id);
	bool pins_migrated = false;
	// transfer ownership
	std::vector<ZOLTAN_ID_TYPE> local_indices(num_ids);
	for(int ii=0; ii < num_ids; ii++) {
		local_indices[ii] = local_ids[ii*num_lid_entries];
	} 

	// DON'T DO MIGRATION (POSSIBLE OUT OF MEMORY IN ARCHER)
	for(int ii=0; ii < num_ids; ii++) {
		// TODO: needs to put data into buf, and index the data in idx array
		//idx_t* partitioning;	// pointer to the partitioning array
		hg->partitioning[global_ids[ii*num_gid_entries]] = dest[ii];
	}
	return;

	// TESTING //
	/////////////
	
	int counts = 0;
	for(int ii=0; ii < num_ids; ii++) {
		ZOLTAN_ID_TYPE global_id = global_ids[ii*num_gid_entries];
		ZOLTAN_ID_TYPE local_id = local_ids[ii*num_lid_entries];
		//for(int jj=0; jj < hg->numMyVertices; jj++) {
		//	if(global_id == hg->vtxGID[jj]) {
		//		counts++;
		//		break;
		//	}
		//}
		//printf("%ld == %ld\n",hg->vtxGID[local_id],global_id);
		//if(hg->process_id == 1) printf("%i: Local %ld global %ld\n",hg->process_id,local_ids[ii*num_lid_entries],global_ids[ii*num_gid_entries]);
		if(hg->vtxGID[local_id] == global_id) counts++;
	}
	printf("%i: VtxGID includes %i out of %i\n",hg->process_id,counts,num_ids);
	/////////////////


	//std::sort(local_indices.begin(), local_indices.end(), [](const int a, const int b) {return a > b; });
	for(int ii=0; ii < num_ids; ii++) {
		// TODO: needs to put data into buf, and index the data in idx array
		//idx_t* partitioning;	// pointer to the partitioning array
		hg->partitioning[global_ids[ii*num_gid_entries]] = dest[ii];
		//ZOLTAN_ID_TYPE local_id = local_ids[ii*num_lid_entries];
		ZOLTAN_ID_TYPE local_id = local_indices[ii];
		//move vertex weight to buffer (represented as int)
		// idx should be indexed based on original object order!
		int current_buf_index = idx[ii];
		current_buf_index += convert_to_bytes((int)(hg->vtx_wgts[local_id]),&(buf[current_buf_index]));
		if(hg->numMyVertices > local_id + 1) {
			//float* vtx_wgts;	/* weights of each of the vertices */
			memmove(&(hg->vtx_wgts[local_id]),&(hg->vtx_wgts[local_id])+1,sizeof(float) * (hg->numMyVertices-local_id-1));
			//ZOLTAN_ID_TYPE *vtxGID;        /* global ID of these vertices */
			memmove(&(hg->vtxGID[local_id]),&(hg->vtxGID[local_id])+1,sizeof(ZOLTAN_ID_TYPE) * (hg->numMyVertices-local_id-1));
		}
		// update local ids
		//std::transform(&(local_indices[ii]),&(local_indices[ii])+(num_ids-ii),&(local_indices[ii]),[local_id](int i) { if(i > local_id) return i-=1; else return i; });	
		std::transform(local_indices.begin()+ii,local_indices.end(),local_indices.begin()+ii,[local_id](int i) { if(i > local_id) return i-=1; else return i; });	
		int pins_for_vertex = hg->vtxedge_ptr[local_id+1] - hg->vtxedge_ptr[local_id];
		if(pins_for_vertex > 0) {
			pins_migrated = true;
			// move data to be migrated
			for(int dd=0; dd < pins_for_vertex; dd++) {
				current_buf_index += convert_to_bytes(hg->pin_GID[hg->vtxedge_ptr[local_id]+dd],&(buf[current_buf_index]));
			}
			if(hg->vtxedge_ptr[local_id]+pins_for_vertex < hg->npins) {
				//ZOLTAN_ID_TYPE *pin_GID;  /* Edges of vertex vtxGID[i] begin at pin_GID[vtxedge_ptr[i]] */ // pin_GID
				memmove(&(hg->pin_GID[hg->vtxedge_ptr[local_id]]),&(hg->pin_GID[hg->vtxedge_ptr[local_id]+pins_for_vertex]),sizeof(ZOLTAN_ID_TYPE) * (hg->npins-(hg->vtxedge_ptr[local_id]+pins_for_vertex)));
			}
			//int *vtxedge_ptr;     /* index into pin_GID array of edge's vertices */ // vtxedge_ptr
			if(hg->numMyVertices > local_id + 1) {
				memmove(&(hg->vtxedge_ptr[local_id]),&(hg->vtxedge_ptr[local_id+1]),sizeof(int) * (hg->numMyVertices-local_id));
				std::transform(&(hg->vtxedge_ptr[local_id]),&(hg->vtxedge_ptr[local_id]) + (hg->numMyVertices-local_id),&(hg->vtxedge_ptr[local_id]),[pins_for_vertex](int i) { return i-=pins_for_vertex; });	
			}
		}
		//int numMyVertices;  /* number of vertices that I own initially */
		hg->numMyVertices--;
		//int npins; /* number of vertices in my hyperedges */
		hg->npins -= pins_for_vertex;
	}
	int pins = 0;
	for(int ii=0; ii < hg->numMyVertices; ii++) {
		pins += hg->vtxedge_ptr[ii+1] - hg->vtxedge_ptr[ii];
	}
	if(pins != hg->npins) printf("%i: npins does not match vtxedge_ptr count of pins\n",hg->process_id);
	if(hg->npins != hg->vtxedge_ptr[hg->numMyVertices]) printf("%i: npins does not match last vtx_ptr value\n",hg->process_id);
	int errors = 0;
	for(int ii=0; ii < num_ids; ii++) {
		ZOLTAN_ID_TYPE global_id = global_ids[ii*num_gid_entries];
		for(int jj=0; jj < hg->numMyVertices; jj++) {
			if(global_id == hg->vtxGID[jj]) {
				errors++;
				break;
			}
		}
	}
	if(errors > 0) printf("%i: Number of vertex errors: %i (out of %i)\n",hg->process_id,errors,num_ids);

	// reallocating memory (pontentially costly!)
	if(num_ids > 0) {
		//hg->vtxGID
		ZOLTAN_ID_TYPE* ptr = (ZOLTAN_ID_TYPE*)realloc(hg->vtxGID,hg->numMyVertices * sizeof(ZOLTAN_ID_TYPE));
		if(ptr != NULL) hg->vtxGID = ptr; 
		else printf("%i: Error realloc in pack function\n",hg->process_id);
		//hg->vtxedge_ptr
		int* ptr2 = (int*)realloc(hg->vtxedge_ptr,(hg->numMyVertices+1) * sizeof(int));
		if(ptr2 != NULL) hg->vtxedge_ptr = ptr2;
		else printf("%i: Error realloc in pack function\n",hg->process_id);
		//hg->vtx_wgts
		float* ptr3 = (float*)realloc(hg->vtx_wgts,hg->numMyVertices * sizeof(float));
		if(ptr3 != NULL) hg->vtx_wgts = ptr3;
		else printf("%i: Error realloc in pack function\n",hg->process_id);
	}
	
	if(pins_migrated) {
		//hg->pin_GID
		ZOLTAN_ID_TYPE* ptr = (ZOLTAN_ID_TYPE*)realloc(hg->pin_GID,hg->npins * sizeof(ZOLTAN_ID_TYPE));
		if(ptr != NULL) hg->pin_GID = ptr;
		else printf("%i: Error realloc in pack function\n",hg->process_id);
	}
	*ierr = ZOLTAN_OK;


	
}
void user_unpack_multi_obj(void *data, int num_gid_entries, int num_ids, ZOLTAN_ID_PTR global_ids, int* sizes, int* idx, char *buf, int *ierr) {
	HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
	PRINTF("%i: user_unpack_multi_obj\n\n",hg->process_id);
	if(num_ids <= 0 ) {
		*ierr = ZOLTAN_OK;
		return;
	}

	// DON'T DO MIGRATION (POSSIBLE OUT OF MEMORY IN ARCHER)
	for(int ii=0; ii < num_ids; ii++) {
		hg->partitioning[global_ids[ii*num_gid_entries]] = hg->process_id;
	}
	return;


	// reallocate memory
	//hg->vtxGID
	ZOLTAN_ID_TYPE* ptr1 = (ZOLTAN_ID_TYPE*)realloc(hg->vtxGID,(hg->numMyVertices+num_ids) * sizeof(ZOLTAN_ID_TYPE));
	if(ptr1 != NULL) hg->vtxGID = ptr1;
	else printf("%i: Error realloc in unpack function\n",hg->process_id);
	//hg->vtxedge_ptr
	int* ptr2 = (int*)realloc(hg->vtxedge_ptr,(hg->numMyVertices+num_ids+1) * sizeof(int));
	if(ptr2 != NULL) hg->vtxedge_ptr = ptr2;
	else printf("%i: Error realloc in unpack function\n",hg->process_id);
	//hg->vtx_wgts
	float* ptr3 = (float*)realloc(hg->vtx_wgts,(hg->numMyVertices+num_ids) * sizeof(float));
	if(ptr3 != NULL) hg->vtx_wgts = ptr3;
	else printf("%i: Error realloc in unpack function\n",hg->process_id);
	// resize hg->pin_GID
	int new_pins = 0;
	for(int ii=0; ii < num_ids; ii++) {
		new_pins += (sizes[ii] - sizeof(int)) / sizeof(ZOLTAN_ID_TYPE);
	}
	if(new_pins > 0) {
		ZOLTAN_ID_TYPE* ptr4 = (ZOLTAN_ID_TYPE*)realloc(hg->pin_GID,(hg->npins+new_pins) * sizeof(ZOLTAN_ID_TYPE));
		if(ptr4 != NULL) hg->pin_GID = ptr4;
		else printf("%i: Error realloc in unpack function\n",hg->process_id);
	}

	

	// add ownership
	int pin_ptr = hg->vtxedge_ptr[hg->numMyVertices];
	
	for(int ii=0; ii < num_ids; ii++) {
		hg->partitioning[global_ids[ii*num_gid_entries]] = hg->process_id;
		//ZOLTAN_ID_TYPE *vtxGID;        /* global ID of these vertices */
		hg->vtxGID[hg->numMyVertices] = global_ids[ii*num_gid_entries];
		//float* vtx_wgts;	/* weights of each of the vertices */
		int current_buf_index = idx[ii];
		int weight;
		current_buf_index += convert_from_bytes(&buf[current_buf_index],&weight);
		hg->vtx_wgts[hg->numMyVertices] = (float)weight;
		hg->vtxedge_ptr[hg->numMyVertices] = pin_ptr;
		//int numMyVertices;  /* number of vertices that I own initially */
		hg->numMyVertices++;
		int pins_for_vertex = (sizes[ii] - sizeof(int)) / sizeof(ZOLTAN_ID_TYPE);
		for(int dd=0; dd < pins_for_vertex; dd++) {
			ZOLTAN_ID_TYPE pin;
			current_buf_index += convert_from_bytes(&buf[current_buf_index],&pin);
			//int npins; /* number of vertices in my hyperedges */
			hg->npins++;
			//ZOLTAN_ID_TYPE *pin_GID;  /* Edges of vertex vtxGID[i] begin at pin_GID[vtxedge_ptr[i]] */ // pin_GID
			hg->pin_GID[pin_ptr] = pin;
			pin_ptr++;
		}
	}
	hg->vtxedge_ptr[hg->numMyVertices] = pin_ptr;
	*ierr = ZOLTAN_OK;

	/////////////
	// testing //
	/////////////
	
	int pins = 0;
	for(int ii=0; ii < hg->numMyVertices; ii++) {
		pins += hg->vtxedge_ptr[ii+1] - hg->vtxedge_ptr[ii];
	}
	if(pins != hg->npins) printf("%i: npins does not match vtxedge_ptr count of pins\n",hg->process_id);
	if(hg->npins != hg->vtxedge_ptr[hg->numMyVertices]) printf("%i: npins does not match last vtx_ptr value\n",hg->process_id);
	int errors = 0;
	for(int ii=0; ii < num_ids; ii++) {
		ZOLTAN_ID_TYPE global_id = global_ids[ii*num_gid_entries];
		bool found = false;
		for(int jj=0; jj < hg->numMyVertices; jj++) {
			if(global_id == hg->vtxGID[jj]) {
				found = true;
				break;
			}
		}
		if(!found) errors++;
	}
	if(errors > 0) printf("%i: Number of vertex errors: %i (out of %i)\n",hg->process_id,errors,num_ids);

	/////////////
	/////////////

	
}
#pragma endregion

class HypergraphPartitioning : public Partitioning {
public:
	
	HypergraphPartitioning(std::vector<Population*>* pops, int population_size) : Partitioning(pops,population_size) {
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
	virtual ~HypergraphPartitioning() {
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
				hg.vtx_wgts[current_local_vertex] = model->null_compute ? 1 : total_vtx_wgts_v[ii];
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
		/* General parameters */
	
		Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
		Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
		Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
		Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
		Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
		Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
		Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* use Zoltan default vertex weights */
		Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */
		Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.005");/* imbalance tolerance */

		// parameters for PHG: http://www.cs.sandia.gov/zoltan/ug_html/ug_alg_phg.html
		
		Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
		Zoltan_Set_Param(zz, "CHECK_HYPERGRAPH", "0");
		Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "MAX");
		Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", "1.0");
		Zoltan_Set_Param(zz, "PHG_CUT_OBJECTIVE", "CONNECTIVITY"); // CONNECTIVITY, HYPEREDGES
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
