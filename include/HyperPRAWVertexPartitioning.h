#ifndef HYPERPRAWVERTEX_PARTITIONING
#define HYPERPRAWVERTEX_PARTITIONING

#include <vector>
#include <algorithm>
#include <set>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>

#include "Simulation.h"
#include "Partitioning.h"
#include "Utils.h"
#include "PRAW.h"

class HyperPRAWVertexPartitioning : public Partitioning {
public:
	
	HyperPRAWVertexPartitioning(std::vector<Population*>* pops, int population_size, char* comm_bandwidth_file, bool max_streams) : Partitioning(pops,population_size) {
		comm_bandwidth_filename = comm_bandwidth_file;
        use_max_streams = max_streams;
	}
	virtual ~HyperPRAWVertexPartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		// fixed parameters for HyperPRAW
        int max_iterations = 15;
        float imbalance_tolerance = 1.2f;
        int num_streams = use_max_streams ? min(96,partitions) : 1; // artificial limit due to memory requirements per stream
        bool local_parallel_update_only = false;
        int sync_batch_size = 1;
        bool input_order_round_robin = true;
        float lambda = 0.8;
        bool save_partitioning_history = true;

        if(partitions <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

        // only active partitions (streams) should bother with it
        // the others will wait for the result and synchronise

        MPI_Comm partitioning_comm;

        // Create MPI communicator group (of only max_number_processes)
        // Only those processes will work towards partitioning
        MPI_Comm_split(MPI_COMM_WORLD,
                        process_id < num_streams ? 1 : 0, // this is the colour, is used to match the processes to the new group
                        process_id,                                         // this is the key, used to generate new rank
                        &partitioning_comm);                               // the subgroup generated
        
        int num_processes;
        MPI_Comm_size(partitioning_comm,&num_processes); 

        bool is_active = process_id < num_streams;

        if(is_active) {
            // PREP DATA STRUCTURES
            // initialise vertex weight values
            int* vtx_wgt = (int*)calloc(model->population_size,sizeof(int));
            for(int ii =0; ii < model->population_size; ii++) {
                vtx_wgt[ii] += 1;
                // Zoltan does consider this balance, so needs to be a fair comparison
                for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
                    vtx_wgt[abs(model->interconnections[ii][jj])] += 1;
                }
            }

            // p2p communication cost estimates from file
            double** comm_cost_matrix = (double**)malloc(sizeof(double*) * partitions);
            for(int ii=0; ii < partitions; ii++) {
                comm_cost_matrix[ii] = (double*)calloc(partitions,sizeof(double));
            }
            // communication cost matrix
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,partitions,false);

            // create hMETIS file (distributed format)
            // each line represents a vertex (neuron) and the ids of the hyperedges (synapsys) that contain it
            // user model->interconnections_size and model->interconnections[][] to build file

            std::string hgraph_file = "model_";
            hgraph_file += "_";
            char str_int[16];
            sprintf(str_int,"%i",num_processes);
            hgraph_file += str_int;
            hgraph_file += "_";
            sprintf(str_int,"%i",process_id);
            hgraph_file += str_int;
            hgraph_file += "_";
            sprintf(str_int,"%i",partitions);
            hgraph_file += str_int;
            hgraph_file += ".hgr";

            int local_stream_size = model->population_size / num_processes + ((model->population_size % num_processes > process_id) ? 1 : 0);
            std::vector<std::vector<int> > hedge_ptr(local_stream_size);
            int counter = 1; // indices start in 1
            for(int ii=0; ii < model->population_size; ii++) {
                // account for presynaptic neuron, as it belongs to the hyperedge too
                int pre_syn_vertex = ii;
                if(pre_syn_vertex % num_processes == process_id) {
                    hedge_ptr[pre_syn_vertex / num_processes].push_back(counter);
                }
                // all postsynaptic connections also belong to the same hyperedge as pre synaptic neuron
                for(int jj=0; jj < model->interconnections_size[pre_syn_vertex]; jj++) {
                    int post_synaptic_vertex = abs(model->interconnections[pre_syn_vertex][jj]);
                    if(post_synaptic_vertex % num_processes == process_id) {
                        hedge_ptr[post_synaptic_vertex / num_processes].push_back(counter);
                    }
                }
                counter++;            
            }

            printf("%i: %li\n",process_id,hedge_ptr.size());

            // create hMETIS (distributed format) file
            PRINTF("%i: Storing model in file %s\n",process_id,hgraph_file.c_str());
            FILE *fp = fopen(hgraph_file.c_str(), "w+");
            
            // write header: NUM_HYPEREDGES NUM_VERTICES
            //  num hyperedges == number of vertices, since each hyperedge represents a presynaptic neuron and all its connecting post synaptic neighbours
            fprintf(fp,"%lu %lu",model->population_size,model->population_size);
            fprintf(fp,"\n");

            for(int ii=0; ii < hedge_ptr.size(); ii++) {
                for(int he=0; he < hedge_ptr[ii].size(); he++) {
                    fprintf(fp,"%i ",hedge_ptr[ii][he]);
                }
                fprintf(fp,"\n");
                
            }
            fclose(fp);
            //////////////////////////////////////////////////////

            // Perform partitioning amongst the streams
            std::string filename = hgraph_file;
            filename += "_prawParallel";
            char experiment_name[filename.length() + 1];
            strcpy(experiment_name, filename.c_str()); 
            int iterations = PRAW::ParallelHDRF(experiment_name, partitioning, partitions, partitioning_comm, comm_cost_matrix, hgraph_file.c_str(), vtx_wgt, max_iterations, imbalance_tolerance,local_parallel_update_only,sync_batch_size,input_order_round_robin,lambda,save_partitioning_history,&previous_activity->at(0));

            // clean up files
            if( remove(hgraph_file.c_str()) != 0 )
                printf( "Error deleting temporary hgraph file %s\n",hgraph_file.c_str() );
            
            // clean memory
            free(vtx_wgt);
            for(int ii=0; ii < partitions; ii++) {
                free(comm_cost_matrix[ii]);
            }
            free(comm_cost_matrix);

        }

        // synchronise partitioning with all processes
        MPI_Barrier(MPI_COMM_WORLD);

        // broadcast results to not participating processes
        MPI_Bcast(partitioning, model->population_size, MPI_LONG, 0,MPI_COMM_WORLD);
	}

private:
    char* comm_bandwidth_filename = NULL;
    bool use_max_streams = false;
};


#endif
