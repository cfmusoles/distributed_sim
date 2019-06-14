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
	
	HyperPRAWVertexPartitioning(std::vector<Population*>* pops, int population_size, char* comm_bandwidth_file, bool useBandwidth) : Partitioning(pops,population_size) {
		comm_bandwidth_filename = comm_bandwidth_file;
        use_bandwidth_file = useBandwidth;
	}
	virtual ~HyperPRAWVertexPartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		int max_iterations = 100;
        float imbalance_tolerance = 1.2f;
        
        if(partitions <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

        // PREP DATA STRUCTURES
        // initialise vertex weight values --> if model->null_compute then uniform; else number of incomming connections
        int* he_wgt = (int*)calloc(model->population_size,sizeof(int));
        // Zoltan does consider this balance, so needs to be a fair comparison
        for(int ii=0; ii < model->population_size; ii++) {
            he_wgt[ii] += 1;
            for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
                he_wgt[abs(model->interconnections[ii][jj])] += 1;
            }
        }

        // p2p communication cost estimates from file
        double** comm_cost_matrix = (double**)malloc(sizeof(double*) * partitions);
        for(int ii=0; ii < partitions; ii++) {
            comm_cost_matrix[ii] = (double*)calloc(partitions,sizeof(double));
        }

        if(use_bandwidth_file)
            PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,partitions,false);
        else
            PRAW::get_comm_cost_matrix_from_bandwidth(NULL,comm_cost_matrix,partitions,false);
        
        // temporary hMETIS file
        std::string hgraph_file = "model_";
        char str_int[16];
        sprintf(str_int,"%i",partitions);
        hgraph_file += str_int;
        hgraph_file += ".hgr";

        // 1 Save current model to hMETIS format file (only one process needs to do this)
        //      Header: num_hyperedges num_vertices
        //      each line after that: a hyperedge which lists the contained vertices
        // 2 Pass on the info to PRAW::ParallelIndependentRestreamingPartitioning (or sequential version) to get partitioning
        // 3 Make sure all partitions have the final partitioning
        // 4 Clean up unused datastructures

        // MASTER NODE creates hMETIS file
        if(process_id == 0) {
            // generate hMETIS file
            // store model in file
            PRINTF("%i: Storing model in file %s\n",process_id,hgraph_file.c_str());
            FILE *fp = fopen(hgraph_file.c_str(), "w+");
            
            // write header: NUM_HYPEREDGES NUM_VERTICES
            //  num hyperedges == number of vertices, since each hyperedge represents a presynaptic neuron and all its connecting post synaptic neighbours
            fprintf(fp,"%lu %lu",model->population_size,model->population_size);
            fprintf(fp,"\n");

            std::vector<std::vector<int> > presynaptics(model->population_size);
            // add presynaptic neuron to hyperedge lead by each post synaptic neuron. See HypergraphPartitioning
            for(int ii=0; ii < model->population_size; ii++) {
                for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
                    presynaptics[abs(model->interconnections[ii][jj])].push_back(ii);
                }
            }

            // write connectivity per neuron id
            for(int ii=0; ii < model->population_size; ii++) {
                // write presynaptic neuron as first neuron
                fprintf(fp,"%i ",ii + 1);
                // write all connections from ii
                for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
                    fprintf(fp,"%i ",abs(model->interconnections[ii][jj]) + 1);
                }
                // write all connections to ii
                /*for(int jj=0; jj < presynaptics[ii].size(); jj++) {
                    fprintf(fp,"%i ",presynaptics[ii][jj] + 1);
                } */
                
                fprintf(fp,"\n");
            }
            fclose(fp);

        } 
        // wait until file is generated
        MPI_Barrier(MPI_COMM_WORLD);

        // define algorithmic parameters
        std::string filename = hgraph_file;
        filename += "_prawParallel";
        char experiment_name[filename.length() + 1]; 
        strcpy(experiment_name, filename.c_str()); 
        PRAW::ParallelHyperedgePartitioning(experiment_name,partitioning, comm_cost_matrix, hgraph_file.c_str(), he_wgt, max_iterations, imbalance_tolerance, false);
        
        if(process_id == 0) {
            if(!use_bandwidth_file && comm_bandwidth_filename != NULL)
                PRAW::get_comm_cost_matrix_from_bandwidth(comm_bandwidth_filename,comm_cost_matrix,partitions,false);
                
            //PRAW::storePartitionStats(filename,partitioning,partitions,model->population_size,&hyperedges,&hedge_ptr,he_wgt,comm_cost_matrix);
            // Use PRAW::getVertexCentricPartitionStatsFromFile instead and record results manually if needed
        }

        free(he_wgt);
        for(int ii=0; ii < partitions; ii++) {
            free(comm_cost_matrix[ii]);
        }
        free(comm_cost_matrix);

        // wait for all before deleting the hgraph file
        MPI_Barrier(MPI_COMM_WORLD);

        if(process_id == 0) {
            // remove graph file
            if( remove(hgraph_file.c_str()) != 0 )
                printf( "Error deleting temporary hgraph file %s\n",hgraph_file.c_str() );
        }

	}

private:
    char* comm_bandwidth_filename = NULL;
    bool use_bandwidth_file = false;
};


#endif