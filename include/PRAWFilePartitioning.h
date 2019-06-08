#ifndef PRAW_PARTITION_PARTITIONING
#define PRAW_PARTITION_PARTITIONING

#include <vector>
#include <algorithm>
#include <set>
#include <cstring>
#include <cstdio>
#include <stdlib.h>

#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"
#include "PRAW.h"

class PRAWFilePartitioning : public Partitioning {
public:
	
	PRAWFilePartitioning(std::vector<Population*>* pops, int population_size, char* comm_bandwidth_file, bool parallel, bool useBandwidth) : Partitioning(pops,population_size) {
		comm_bandwidth_filename = comm_bandwidth_file;
        isParallel = parallel;
        use_bandwidth_file = useBandwidth;
	}
	virtual ~PRAWFilePartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		int max_iterations = 100;
        float imbalance_tolerance = 1.1f;
        
        if(partitions <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

        // Assign unique vertices to partitions
        for(int ii=0; ii < model->population_size; ii++) {
            partitioning[ii] = (double)rand() / (double)RAND_MAX * partitions;
            if(partitioning[ii] == partitions) partitioning[ii] -= 1;
        }
        // 1 - Convert graph into hMETIS format
        //      list of hyperedges (each with a list of vertex id)
        //      list of vertex pointers (each with a list of the hyperedges a vertex belongs to)
        // 2 - Pass data to PRAW to partition
        // 3 - return partitioning
        
        /*if(model->hypergraph_file == NULL) {
            PRINTF("%i: hypergraph file not set in Model object. Random partition.",process_id);
            return;
        }

        std::vector<std::vector<int> > hyperedges(model->population_size);
        std::vector<std::vector<int> > hedge_ptr(model->population_size);
        
        PRAW::load_hypergraph_from_file(model->hypergraph_file,&hyperedges,&hedge_ptr);*/
        

        // initialise vertex weight values --> if model->null_compute then uniform; else number of incomming connections
        int* vtx_wgt = (int*)calloc(model->population_size,sizeof(int));
        for(int ii =0; ii < model->population_size; ii++) {
            vtx_wgt[ii] = 1;
        }

        /*if(!model->null_compute) {
            // weight of computation associated to incoming synapses
            // assume all to all connectivity within a hyperedge (needs to match how model->interconnectivity is constructed)
            for(int ii=0; ii < model->population_size; ii++) {
                for(int jj=0; jj < hedge_ptr[ii].size(); jj++) {
                    int he_id = hedge_ptr[ii][jj];
                    for(int kk=0; kk < hyperedges[he_id].size(); kk++) {
                        int to = hyperedges[he_id][kk];
                        vtx_wgt[to] += hyperedges[he_id].size()-1;
                    }
                }
            }
        }*/

        // p2p communication cost estimates from file
        double** comm_cost_matrix = (double**)malloc(sizeof(double*) * partitions);
        for(int ii=0; ii < partitions; ii++) {
            comm_cost_matrix[ii] = (double*)calloc(partitions,sizeof(double));
        }

        if(use_bandwidth_file)
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,partitions,false);
        else
            PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,partitions,false);
        
        std::string filename = model->hypergraph_file;
        if(isParallel) {
            filename += "_prawParallel";
            // define parameters
            float ta_refine = 0.95f;
            int stopping_condition = 2;
            char* experiment_name;
            PRAW::ParallelIndependentRestreamingPartitioning(experiment_name,partitioning,comm_cost_matrix, model->hypergraph_file,vtx_wgt,max_iterations, imbalance_tolerance,ta_refine,true,stopping_condition,false);
        } else {
            filename += "_prawSequential";
            // define parameters
            float ta_refine = 0.95f;
            int stopping_condition = 2;
            char* experiment_name;
            PRAW::SequentialStreamingPartitioning(experiment_name,partitioning,partitions,comm_cost_matrix, model->hypergraph_file,vtx_wgt,max_iterations, imbalance_tolerance,ta_refine,true,stopping_condition, false);
        }

        if(process_id == 0) {
            if(!use_bandwidth_file && comm_bandwidth_filename != NULL)
                PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,partitions,false);
                
            //PRAW::storePartitionStats(filename,partitioning,partitions,model->population_size,&hyperedges,&hedge_ptr,vtx_wgt,comm_cost_matrix);
            // Use PRAW::getVertexCentricPartitionStatsFromFile instead and record results manually if needed
        }

        free(vtx_wgt);
        for(int ii=0; ii < partitions; ii++) {
            free(comm_cost_matrix[ii]);
        }
        free(comm_cost_matrix);

	}

private:
    char* comm_bandwidth_filename = NULL;
    bool isParallel = false;
    bool use_bandwidth_file = false;
};


#endif
