#ifndef HYPERPRAWEDGE_PARTITIONING
#define HYPERPRAWEDGE_PARTITIONING

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

class HyperPRAWEdgePartitioning : public Partitioning {
public:
	
	HyperPRAWEdgePartitioning(std::vector<Population*>* pops, int population_size, char* comm_bandwidth_file, bool parallel, bool useBandwidth) : Partitioning(pops,population_size) {
		comm_bandwidth_filename = comm_bandwidth_file;
        isParallel = parallel;
        use_bandwidth_file = useBandwidth;
	}
	virtual ~HyperPRAWEdgePartitioning() {}
	
	virtual void perform_partitioning(Model* model, int partitions, int process_id, std::vector<int>* previous_activity) {
		int max_iterations = 100;
        float imbalance_tolerance = 1.1f;
        
        if(partitions <= 1) {
			PRINTF("Partitioning not required\n");
			return;
		}

        // PREP DATA STRUCTURES
        // initialise vertex weight values
        int* vtx_wgt = (int*)calloc(model->population_size,sizeof(int));
        for(int ii =0; ii < model->population_size; ii++) {
            vtx_wgt[ii] = 1;
        }
        // Zoltan does consider this balance, so needs to be a fair comparison
        for(int ii=0; ii < model->population_size; ii++) {
            for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
                vtx_wgt[abs(model->interconnections[ii][jj])] += 1;
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
        
        
        // 1 Save current model to hMETIS format file (only one process needs to do this)
        //      Header: num_hyperedges num_vertices
        //      each line after that: a hyperedge which lists the contained vertices
        // 2 Pass on the info to PRAW::ParallelIndependentRestreamingPartitioning (or sequential version) to get partitioning
        // 3 Make sure all partitions have the final partitioning
        // 4 Clean up unused datastructures

        /////////////////////////////////////////
        //for PRAW::ParallelHyperedgePartitioning
        /////////////////////////////////////////
        // temporary hMETIS file
        /*std::string hgraph_file = "model_";
        char str_int[16];
        sprintf(str_int,"%i",partitions);
        hgraph_file += str_int;
        hgraph_file += ".hgr";
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

            // write connectivity per neuron id
            for(int ii=0; ii < model->population_size; ii++) {
                // write presynaptic neuron as first neuron
                fprintf(fp,"%i ",ii + 1);
                // write all connections from ii
                for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
                    fprintf(fp,"%i ",abs(model->interconnections[ii][jj]) + 1);
                }
                
                fprintf(fp,"\n");
            }
            fclose(fp);

        }
        ////////////////////////////////////////
         */


        /////////////////////////////////////////
        //for PRAW::ParallelHDRF
        /////////////////////////////////////////
        std::string hgraph_file = "model_";
        hgraph_file += "_";
        char str_int[16];
        sprintf(str_int,"%i",partitions);
        hgraph_file += str_int;
        hgraph_file += "_";
        sprintf(str_int,"%i",process_id);
        hgraph_file += str_int;
        hgraph_file += ".hgr";
        // create hMETIS file (distributed format)
        // each line represents a vertex
        // the contents of the line are the hyperedge ids that vertex belongs to
        // generate intermediate datastructure with list of hedge ids

        std::vector<std::vector<int> > hedge_ptr(ceil((float)model->population_size / partitions));
        // read reminder of file (one line per hyperedge)
        int counter = 1;
        for(int ii=0; ii < model->population_size; ii++) {
            // account for presynaptic neuron, as it belongs to the hyperedge too
            int pre_syn_vertex = ii;
            if(pre_syn_vertex % partitions == process_id) {
                hedge_ptr[pre_syn_vertex / partitions].push_back(counter);
            }
            // all postsynaptic connections
            for(int jj=0; jj < model->interconnections_size[ii]; jj++) {
                int vertex_id = abs(model->interconnections[ii][jj]);
                if(vertex_id % partitions == process_id) {
                    hedge_ptr[vertex_id / partitions].push_back(counter);
                }
            }
            counter++;            
        }

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



        MPI_Barrier(MPI_COMM_WORLD); 

        // define algorithmic parameters
        float ta_refine = 0.95f;
        int stopping_condition = 2;
        std::string filename = hgraph_file;
        if(isParallel) {
            filename += "_prawParallel";
            char experiment_name[filename.length() + 1]; 
            strcpy(experiment_name, filename.c_str()); 
            //PRAW::ParallelHyperedgePartitioning(experiment_name,partitioning,comm_cost_matrix, hgraph_file.c_str(),vtx_wgt,max_iterations, imbalance_tolerance,ta_refine,true,stopping_condition,false);
            // alternative method (based on Alistairh streaming, only one vertex with all its hyperedges at a time)
            //PRAW::ParallelHyperedgePartitioning_he_stream(experiment_name,partitioning,comm_cost_matrix, model->population_size,&he_stream,vtx_wgt,max_iterations, imbalance_tolerance,false);
            PRAW::ParallelHDRF(experiment_name,partitioning, comm_cost_matrix, hgraph_file, vtx_wgt, max_iterations, imbalance_tolerance, true);
        } else {
            filename += "_prawSequential";
            char experiment_name[filename.length() + 1]; 
            strcpy(experiment_name, filename.c_str()); 
            PRAW::SequentialHyperedgePartitioning(experiment_name,partitioning,partitions,comm_cost_matrix, hgraph_file.c_str(),vtx_wgt,max_iterations, imbalance_tolerance,ta_refine,true,stopping_condition, false);
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

        // wait for all before deleting the hgraph file
        MPI_Barrier(MPI_COMM_WORLD);

        // remove hMETIS files
        /*if(process_id == 0) {
            // remove graph file
            if( remove(hgraph_file.c_str()) != 0 )
                printf( "Error deleting temporary hgraph file %s\n",hgraph_file.c_str() );
        }*/
        if( remove(hgraph_file.c_str()) != 0 )
            printf( "Error deleting temporary hgraph file %s\n",hgraph_file.c_str() );

	}

private:
    char* comm_bandwidth_filename = NULL;
    bool isParallel = false;
    bool use_bandwidth_file = false;
};


#endif
