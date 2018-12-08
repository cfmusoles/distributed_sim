//#define API_PROFILING					// Scalasca User API profiling
//#define RECORD_STATISTICS				// Store network activity and stats in files
//#define RECORD_ACTIVITY					// store neuron activity in file
//#define RECORD_PROCESSES_ACTIVITY		// store processes timings (one file per process)
//#define RECORD_PROCESS_TRACE			// store processes trace information (timing for compute, comm, etc.)
//#define MEASURE_IDLE_TIME				// Separate comm into process idle (wait for others) and sync time
//#define ADVANCED_COMM_STATS				// store, for each communication, size of message
#define ADVANCED_COMM_STATS_MATRIX_ONLY	// only store process-to-process total messaging (not individual message size) 
//#define VERBOSE						// display debugging information (build time)
#define PARTITION_CONNECTIVITY_GRAPH	// store partition connectivity graph
					

// TODO
// Communicator.neuron_connection_to should not store lists for non local neurons
// check any memory footprint that scales with the number of processes

#include <utmpx.h>
#include <metis.h>
#include "parmetis.h"
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <numeric>
#include <typeinfo>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <json/json.h>

#include "Neuron.h"
#include "Synapse.h"
#include "CurrentInjector.h"
#include "PoissonInjector.h"
#include "Simulation.h"
#include "Utils.h"
#include "NeuronParams.h"
#include "AllGatherCommunicator.h"
#include "DirectCommunicator.h"
#include "SubscriberCommunicator.h"
#include "PEXCommunicator.h"
#include "RMACommunicator.h"
#include "NBXCommunicator.h"
#include "PCXCommunicator.h"
#include "RSXCommunicator.h"
#include "GraphPartitionPartitioning.h"
#include "ParallelPartitionPartitioning.h"
#include "HypergraphPartitioning.h"
#include "HIERPartitioning.h"
#include "RandomPartitioning.h"
#include "RandomBalancedPartitioning.h"
#include "CustomPartitioning.h"
#include "PRAWFilePartitioning.h"
#include "ZoltanFilePartitioning.h"
#include "RoundRobinPartitioning.h"
#include "Connectivity.h"
#include "RandomConnectivity.h"
#include "StaticConnectivity.h"


#ifdef API_PROFILING
#include <scorep/SCOREP_User.h>
#endif


const int MPI_PROCESS_WAIT_AHEAD = 0.0f;		// time each MPI process waits once spawn before processing or start counting
const int MASTER_NODE = 0;

bool MODEL_LOADED_FROM_FILE = false;	// internal flag to indicate current model has been loaded from file

bool is_MPI_error(int error_code) {
	if(error_code == MPI_SUCCESS) {
		return false;
	} else {
		char error_string[BUFSIZ];
		int error_length;
		MPI_Error_string(error_code, error_string, &error_length);
		printf("********Error code: %s\n",error_string);
		return true;
	}
}

void exit_error() {
	MPI_Finalize();
	exit(-1);
}


struct IdxCompare
{
    const std::vector<Synapse>& target;

    IdxCompare(const std::vector<Synapse>& target): target(target) {}

    bool operator()(int a, int b) const { return target[a] < target[b]; }
};

// Load connectivity from file
// FORMAT
// HEADER: TOTAL P1 P2... Number of vertices in total, followed by a value per population (sign indicates Exc,+, or Inh,-)
// EACH LINE --> a vertex
// 		VERTEX_WEIGHT PodID E1 E2.... where E1 is the vertex id of the vertex it is connected to
// TODO: Injectors must be stored as well (injector connectivity)
// ERROR with fruit fly loaded model, pop_counter[pop_id] goes beyond the number of cells in a population for the fruit fly loaded model
		
void load_model_from_python_file(std::string filename, NeuronParams* c_params, std::vector<Population*>* pops, std::vector<Connectivity*>* conns) {
	MODEL_LOADED_FROM_FILE = true;
	
	std::ifstream istream(filename.c_str());
	
	if(!istream) {
		printf("Error while opening model file %s\n",filename.c_str());
		return;
	}
	
	std::string line;
	// process header
	std::getline(istream,line);
	std::istringstream buf(line);
    std::istream_iterator<int> beg(buf), end;
	std::vector<int> tokens(beg, end);
	int total_vertices = tokens[0];
	for(int ii=1; ii < tokens.size(); ii++) {
		if(tokens[ii] == 0) 
			continue;
		Population_type pop_type = tokens[ii] > 0 ? EXC : INH;
		char buf[BUFSIZ];
		snprintf(buf, sizeof(buf), "Pop%d", ii);
		Population* p = new Population(buf,abs(tokens[ii]),pop_type,c_params);
		pops->push_back(p);
		StaticConnectivity* c = new StaticConnectivity(p);
		c->connections.resize(abs(tokens[ii])); //total_vertices;
		conns->push_back(c);
	}

	PRINTF("Loaded from file: Populations: %lu; connections %lu:\n",pops->size(),conns->size());
	// read reminder of file (one line per vertex)
	int counter = 0;
	int* pop_counter = (int*)calloc(pops->size(),sizeof(int));
	while(std::getline(istream,line)) {
		//if(rand01() < 0.00) continue;
		std::istringstream buf(line);
		std::istream_iterator<int> beg(buf), end;
		std::vector<int> tokens(beg, end);
		int pop_id = tokens[0];
		Connectivity* c = (*conns)[pop_id];
		(*pops)[pop_id]->neuron_ids.push_back(counter);
		// pop_counter[pop_id] goes beyond the number of cells in a population for the fruit fly loaded model
		for(int ii=1; ii < tokens.size(); ii++) { // skip vertex weight for now
			c->connections[pop_counter[pop_id]].push_back(tokens[ii]);
		}
		counter++;
		pop_counter[pop_id] += 1;
	}
	
	free(pop_counter);
	istream.close();
	
}

// load model from hMETIS formatted file
// first line indicates number of hyperedges and number of vertices
// thereafter, list of vertex id belonging to a hyperedge (per line)
void load_model_from_hmetis_file(char* filename, NeuronParams* c_params, std::vector<Population*>* pops, std::vector<Connectivity*>* conns) {
	MODEL_LOADED_FROM_FILE = true;
	
	std::ifstream istream(filename);
	
	if(!istream) {
		printf("Error while opening hMETIS file %s\n",filename);
		return;
	}

    std::string line;
	// process header
	std::getline(istream,line);
	std::istringstream buf(line);
    std::istream_iterator<int> beg(buf), end;
	std::vector<int> tokens(beg, end);
	int total_vertices = tokens[1];
	int total_hyperedges = tokens[0];

	// all vertices in one population
    Population* p = new Population("Pop1",total_vertices,EXC,c_params);
	pops->push_back(p);
	StaticConnectivity* c = new StaticConnectivity(p);
	c->connections.resize(total_vertices);
	conns->push_back(c);
	for(int ii=0; ii < total_vertices; ii ++) {
		(*pops)[0]->neuron_ids.push_back(ii);
	}
    
	PRINTF("Loaded from hmetis file: Vertices: %i; hyperedges %i:\n",total_vertices,total_hyperedges);
	
    // read reminder of file (one line per hyperedge)
	while(std::getline(istream,line)) {
		std::istringstream buf(line);
		std::istream_iterator<int> beg(buf), end;
		std::vector<int> tokens(beg, end);
		// in each hyperedge, all to all connectivity
		for(int ii=0; ii < tokens.size(); ii++) {
            int vertex_id = tokens[ii]-1;
			for(int jj=ii+1; jj < tokens.size(); jj++) {
				//if(std::find((*conns)[0]->connections[vertex_id].begin(), 
				//		(*conns)[0]->connections[vertex_id].end(),
				//		 tokens[jj]-1) == (*conns)[0]->connections[vertex_id].end()) {
					(*conns)[0]->connections[vertex_id].push_back(tokens[jj]-1);
				//}
			}
			
		}
	}
	istream.close();
}

/* GENERAL WEAKNESSES:
 * Memory requirements per process: 
 * 		an array of Neuron of size number of neurons
 * 		an array of vector<Synapse> of size number of neurons
 * 		partitioning arrays: connections (number of neurons X number of synapses), partitioning (number of neurons)
 */


/* TODO
 * Use CSR distributed file format from ParMETIS to define graph distributed (less mem requirements)
 * Separate simulator from model
 */

int main(int argc, char** argv) {

	if(argc <= 3) {
		printf("Name of the experiment required.\nUsage: distSim -n [name] -c [comm_pattern] -p [partitioning] -s [seed] -w [neuron_activity_file] -N (for null compute) -h [hMETIS_GRAPH] -b [comm_cost_matrix_file]");
		return 0;
	}

	char* sim_name;// = argv[1];
	char* comm_pattern;// = argv[2];
	char* partitioning;// = argv[3];
	char* comm_bandwidth_matrix_file = NULL;
	int r_seed;
	char* w_file;
	bool use_seed = false;
	bool use_weight_file = false;
	float n_scale = 1.0;
	float k_scale = 1.0;
	int t_end = 100;
	char* model_selected = NULL;
	int node_size = 24;

	Model model;
	model.store_in_file = false;
	model.null_compute = false;
	model.hypergraph_file = NULL;

	// getting command line parameters
	extern char *optarg;
	extern int optind, opterr, optopt;
	int c;
	while( (c = getopt(argc,argv,"n:c:p:s:w:Nh:b:k:f:t:m:i:")) != -1 ) {
		switch(c) {
			case 'n': // test name
				sim_name = optarg;
				break;
			case 'c': // communication pattern
				comm_pattern = optarg;
				break;
			case 'p': // partitioning
				partitioning = optarg;
				break;
			case 's': // random seed
				r_seed = atoi(optarg);
				use_seed = true;
				break;
			case 'w': // neuron activity file (weights)
				w_file = optarg;
				use_weight_file = true;
				break;
			case 'N':
				model.null_compute = true;
				break;
			case 'h':
				model.hypergraph_file = optarg;
				break;
			case 'b':
				comm_bandwidth_matrix_file = optarg;
				break;
			case 'k': // fraction of synapses
				k_scale = 0.001 * atoi(optarg);
				break;
			case 'f': // fraction of neurons
				n_scale = 0.001 * atoi(optarg);
				break;
			case 't': // simulated time in ms
				t_end = atoi(optarg);
				break;
			case 'm': // model selection
				model_selected = optarg;
				break;
			case 'i': // processes per node
				node_size = atoi(optarg);
				break;
		}
	}

	if(model.hypergraph_file == NULL && model_selected == NULL) {
		printf("Error, model has not been selected. Use -m parameter to select model ('-m cm' or '-m mcv'\n");
		MPI_Finalize();
		return 0;
	}
	PRINTF("Model selected %s, Neuron scale: %f, Synaptic scale: %f, Simulated time: %i\n",model_selected,n_scale,k_scale,t_end);

	// gathering comm patterns and partitioning methods
	// user can supply multiple via : syntax (example: -c nbx:pex)
	// when supplying multiple methods, simulations are run once per permutation of partitioning / comm pattern
	char * pch;
	pch = strtok (partitioning,":");
	std::vector<std::string> partitioning_methods;
	while (pch != NULL)
	{
		partitioning_methods.push_back(std::string(pch));
		pch = strtok (NULL, ":");
	}
	pch = strtok (comm_pattern,":");
	std::vector<std::string> comm_pattern_methods;
	while (pch != NULL)
	{
		comm_pattern_methods.push_back(std::string(pch));
		pch = strtok (NULL, ":");
	}

	
	MPI_Init(&argc,&argv);
	
	// size of the communicator (number of processes)
	int num_processes;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
	
	// process rank
	int process_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
		
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    PRINTF("Running from processor %s (%i), rank %d"
           " out of %d processes\n",
           processor_name, sched_getcpu(),process_id, num_processes);
	
	
#ifdef API_PROFILING
	SCOREP_RECORDING_OFF();
#endif

	// introduce artificial wait of 60s (ShARC workaround fix suggested by Anthony)
    float the_time = MPI_Wtime();
    do {
            //wait
    } while(MPI_Wtime() - the_time <= MPI_PROCESS_WAIT_AHEAD);

	// broadcast random seed so all partitions work with the same random generator
	int random_seed;
	if(process_id == MASTER_NODE) {
		printf("Running experiment: %s\n\n",sim_name);
		random_seed = time(NULL);
	}
	
	if(use_seed) {
		random_seed = r_seed;
	} else {
		MPI_Bcast(&random_seed, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
	}
	PRINTF("Process %i seed %i\n",process_id,random_seed);
	srand(random_seed);
	
	// cap node_size to num_processes maximum
	if(node_size > num_processes) node_size = num_processes;

	float dt = 0.1f; 					// timestep (ms)
	//int t_end = 350; 					// simulation time (ms)
	
	std::vector<Population*> pops;
	std::vector<Connectivity*> conns;
	std::vector<Injector*> injections;
	std::vector<NeuronParams*> cell_params;
	idx_t* custom_partitioning = NULL;
	
	// USER INPUT //
	// NETWORK CREATION //

	//float n_scale = 0.20f; // scale number of neurons
	//float k_scale = 1.0f;	// scale power of synapses
	
	if(model.hypergraph_file != NULL) {
		// FROM FILE: hMETIS format
		NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
		// no need to set cell params (required only if !model.null_compute)
		cell_params.push_back(c_params);
		c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
		c_params->v_start_dev = 5.0f;		// std for start voltage (mV) was 5
		c_params->tau_m = 10.0f;			// membrane time constant (ms) 
		c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
		c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
		c_params->v_rest = -45.0f;			// resting potential (mV)
		c_params->v_thres = -50.0f;			// spike threshold (mV) 
		c_params->v_reset = -65.0f;			// voltage after spike (mV) 
		c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
		c_params->i_offset = 0.95f;			// constant current input (nA)
		c_params->c_m = 0.250f;				// cm (nF)
		c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
		c_params->w_exc_dev = 0.0088f;			// nA (PyNN default is nA) was 0.0088
		c_params->w_inh_g = -4.0f;			// g is a multiplier factor with respect to W_EXC
		c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
		c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
		c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms)
		c_params->exc_delay_dev = 0.7f;		// delay on exc synaptic propagation (ms) was 0.7
		c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
		c_params->inh_delay_dev = 0.4f;		// delay on inh synaptic propagation (ms) was 0.4
		
		load_model_from_hmetis_file(model.hypergraph_file,c_params,&pops,&conns);
	} else if(strcmp(model_selected,"cm") == 0) {
		// USED FOR FRONTIERS IN NEUROINFORMATICS PAPER //
		// CORTICAL MICROCIRCUIT without injectors. To generate activity:
		// v_rest > v_threshold
		// constant input (i_offset)
		float conn_threshold = 0.0f;		// cap connectivity probability (anything below this value will be set to 0)
		
		NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
		cell_params.push_back(c_params);
		// Values from Potjans and Diesmann 2014 paper
		c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
		c_params->v_start_dev = 5.0f;		// std for start voltage (mV) was 5
		c_params->tau_m = 10.0f;			// membrane time constant (ms) 
		c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
		c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
		c_params->v_rest = -45.0f;			// resting potential (mV)
		c_params->v_thres = -50.0f;			// spike threshold (mV) 
		c_params->v_reset = -65.0f;			// voltage after spike (mV) 
		c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
		c_params->i_offset = 0.95f;			// constant current input (nA)
		c_params->c_m = 0.250f;				// cm (nF)
		c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
		c_params->w_exc_dev = 0.0088f;			// nA (PyNN default is nA) was 0.0088
		c_params->w_inh_g = -4.0f;			// g is a multiplier factor with respect to W_EXC
		c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
		c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
		c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms)
		c_params->exc_delay_dev = 0.7f;		// delay on exc synaptic propagation (ms) was 0.7
		c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
		c_params->inh_delay_dev = 0.4f;		// delay on inh synaptic propagation (ms) was 0.4
			
		// populations
		int n_full[8] = {20683,5834,21915, 5479,4850,1065,14395,2948};
		
		Population* l23e = new Population("l2_3e",n_full[0] * n_scale,EXC,c_params); pops.push_back(l23e);
		Population* l23i = new Population("l2_3i",n_full[1] * n_scale,INH,c_params); pops.push_back(l23i);
		Population* l4e = new Population("l4e",n_full[2] * n_scale,EXC,c_params); pops.push_back(l4e);
		Population* l4i = new Population("l4i",n_full[3] * n_scale,INH,c_params); pops.push_back(l4i);
		Population* l5e = new Population("l5e",n_full[4] * n_scale,EXC,c_params); pops.push_back(l5e);
		Population* l5i = new Population("l5i",n_full[5] * n_scale,INH,c_params); pops.push_back(l5i);
		Population* l6e = new Population("l6e",n_full[6] * n_scale,EXC,c_params); pops.push_back(l6e);
		Population* l6i = new Population("l6i",n_full[7] * n_scale,INH,c_params); pops.push_back(l6i);
		
		// connections ([target][source])
		float conn_probs[8][8] = {{0.1009f,  0.1689f, 0.0437f, 0.0818f, 0.0323f, 0.f, 0.0076f, 0.f},
								{0.1346f,   0.1371f, 0.0316f, 0.0515f, 0.0755f, 0.f,     0.0042f, 0.f},
								{0.0077f,   0.0059f, 0.0497f, 0.135f,  0.0067f, 0.0003f, 0.0453f, 0.f},
								{0.0691f,   0.0029f, 0.0794f, 0.1597f, 0.0033f, 0.f,     0.1057f, 0.f},
								{0.1004f,   0.0622f, 0.0505f, 0.0057f, 0.0831f, 0.3726f, 0.0204f, 0.f},
								{0.0548f,   0.0269f, 0.0257f, 0.0022f, 0.06f,   0.3158f, 0.0086f, 0.f},
								{0.0156f,   0.0066f, 0.0211f, 0.0166f, 0.0572f, 0.0197f, 0.0396f, 0.2252f},
								{0.0364f,   0.001f,  0.0034f, 0.0005f, 0.0277f, 0.008f,  0.0658f, 0.1443f}};
		
		// cap connection probabilities
		for(int ii=0; ii < 8; ii++) {
			for(int jj=0; jj < 8; jj++) {
				if(conn_probs[ii][jj] < conn_threshold) conn_probs[ii][jj] = 0;
			}
		}
			
		// in-degree of connectivity
		float k_full[8][8] = {0};
		for(int ii=0; ii < pops.size(); ii++) {
			for(int jj=0; jj < pops.size(); jj++) {
				k_full[ii][jj] = roundf(log(1.f-conn_probs[ii][jj])/(log((double)((n_full[ii] * n_full[jj] - 1)) / (double)((n_full[ii] * n_full[jj]))))) / n_full[ii];
			}
		}
		// build connections
		for(int ii = 0; ii < pops.size(); ii++) {
			for(int jj=0; jj < pops.size(); jj++) {
				Connectivity* c = new RandomConnectivity(pops[ii],pops[jj],k_scale * k_full[jj][ii] * pops[jj]->num_neurons / (float)(pops[ii]->num_neurons * pops[jj]->num_neurons)); //conn_probs[jj][ii] * k_scale);
				conns.push_back(c);
			}
		}
	} else if(strcmp(model_selected,"mvc") == 0) {
		// USED FOR FRONTIERS IN NEUROINFORMATICS PAPER //
		// MULTI AREA MODEL --> 32 areas, each one similar connectivity CM
		// Schmidt 2018 https://github.com/INM-6/multi-area-model
		// Neurons: ~4.13 million * n_scale
		// Synapses: ~24.2 billion * n_scale^2 * k_scale
		// Model only uses type I (intra area) and III (coricocortical) synapses
		// No external input 
			// constant input (i_offset)
		// Cell params from CM
		NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
		cell_params.push_back(c_params);
		// Values from Potjans and Diesmann 2014 paper
		c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
		c_params->v_start_dev = 5.0f;		// std for start voltage (mV) was 5
		c_params->tau_m = 10.0f;			// membrane time constant (ms) 
		c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
		c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
		c_params->v_rest = -65.0f;			// resting potential (mV)
		c_params->v_thres = -50.0f;			// spike threshold (mV) 
		c_params->v_reset = -65.0f;			// voltage after spike (mV) 
		c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
		c_params->i_offset = 0.38f;			// constant current input (nA)
		c_params->c_m = 0.250f;				// cm (nF)
		c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
		c_params->w_exc_dev = 0.0088f;			// nA (PyNN default is nA) was 0.0088
		c_params->w_inh_g = -15.0f;			// g is a multiplier factor with respect to W_EXC was -4
		c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
		c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
		c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms)
		c_params->exc_delay_dev = 0.7f;		// delay on exc synaptic propagation (ms) was 0.7
		c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
		c_params->inh_delay_dev = 0.4f;		// delay on inh synaptic propagation (ms) was 0.4

		// json examples and doc https://en.wikibooks.org/wiki/JsonCpp
		std::ifstream ifs("multiareaModel.json");
		Json::Reader reader;
		Json::Value root;
		reader.parse(ifs, root); // reader can also read strings

		// population sizes in root["neuron_numbers"][AREA][POPULATION]
		// population should be modulated by n_scale
		std::map<std::string,Population*> mappedPops;
		const Json::Value& areas = root["neuron_numbers"]; // array
		std::vector<std::string> areasKeys = areas.getMemberNames();
		int ns = 0;
		for(size_t a=0; a < areasKeys.size(); a++) {
			const std::string& areaKey = areasKeys[a];
			const Json::Value& populations = areas[areaKey];
			std::vector<std::string> populationsKeys = populations.getMemberNames();
			for(size_t p=0; p < populationsKeys.size(); p++) {
				const std::string& popKey = populationsKeys[p];
				float popSize = populations[popKey].asFloat();
				// pop name --> areaKey + popKey
				// pop size --> (int)popSize.asFloat() * n_scale
				// pop type if pop name contains e = EXC, else INH
				// ignore 'total' values
				std::string e_string = "total";
				// discard "total" values and empty population
				if(popKey.find(e_string) != std::string::npos || popSize <= 0) continue;
				e_string = "E";
				std::string pName = areaKey+popKey;
				Population* pop = new Population(pName.c_str(),
											popSize * n_scale,
											popKey.find(e_string) != std::string::npos ? EXC : INH,
											c_params);
				pops.push_back(pop);
				mappedPops[pName] = pop;
				ns += popSize * n_scale;
			}
		}

		// connectivity (number of synapses) in root["synapses"][target_area][target_pop][source_area][source_pop]
		// connectivity should be modulated based on n_scale and k_scale
		const Json::Value& targetAreas = root["synapses"]; // array
		std::vector<std::string> targetAreaKeys = targetAreas.getMemberNames();
		long int syns = 0;
		for(size_t ta=0; ta < targetAreaKeys.size(); ta++) {
			const std::string& taKey = targetAreaKeys[ta];
			const Json::Value& targetPops = targetAreas[taKey];
			if(targetPops.type() != Json::objectValue) continue;
			std::vector<std::string> targetPopsKeys = targetPops.getMemberNames();
			for(size_t tp=0; tp < targetPopsKeys.size(); tp++) {
				const std::string& tpKey = targetPopsKeys[tp];
				const Json::Value& sourceAreas = targetPops[tpKey];
				if(sourceAreas.type() != Json::objectValue) continue;
				std::vector<std::string> sourceAreasKeys = sourceAreas.getMemberNames();
				for(size_t sa=0; sa < sourceAreasKeys.size(); sa++) {
					const std::string& saKey = sourceAreasKeys[sa];
					const Json::Value& sourcePops = sourceAreas[saKey];
					if(sourcePops.type() != Json::objectValue) continue;
					std::vector<std::string> sourcePopsKeys = sourcePops.getMemberNames();
					for(size_t sp=0; sp < sourcePopsKeys.size(); sp++) {
						const std::string& spKey = sourcePopsKeys[sp];
						if(sourcePops[spKey].type() != Json::realValue) continue;
						float number_synapses = sourcePops[spKey].asFloat();
						//printf("[%s][%s][%s][%s] = %f\n",taKey.c_str(),tpKey.c_str(),saKey.c_str(),spKey.c_str(),sourcePops[spKey].asFloat());
						// connectivity: 
							// from taKey + tpKey --> saKey + spKey
							// number synapses --> sourcePops[spKey].asFloat()
						std::string sourcePop = saKey + spKey;
						std::string targetPop = taKey + tpKey;
						if(mappedPops.count(sourcePop) <= 0 || mappedPops.count(targetPop) <= 0) continue; // population was not defined in the population list
						int source_neurons = mappedPops[sourcePop]->num_neurons;
						int target_neurons = mappedPops[targetPop]->num_neurons;
						float ratio = number_synapses / (float)(source_neurons / n_scale * target_neurons / n_scale);
						float p_conn = k_scale * ratio;
						if(p_conn <= 0) {
							continue;
						}
						//printf("%f (%i,%i) --> (%f) final %f\n",number_synapses,source_neurons,target_neurons ,ratio,p_conn);
						// create connection between populations
						Connectivity* c = new RandomConnectivity(
											mappedPops[sourcePop],
											mappedPops[targetPop],
											p_conn
						); 
						conns.push_back(c);
						syns += (long int)((long int)(source_neurons * target_neurons) * p_conn);
					}
				}
			}
		}
		
		PRINTF("From model --> Number of neurons: %i. Approx synapses: %li\n",ns,syns);
	} else {
		PRINTF("%i: Custom model being used\n",process_id);

	}
	
	 
	
	/*
	// Cortical microcircuit, from (Potjans and Diesman, 2014) and (Urgese, 2016)
	//
	// Parameters taken from https://github.com/NeuralEnsemble/PyNN/tree/master/examples/Potjans2014
	// Parameters taken from https://github.com/NeuralEnsemble/PyNN/tree/master/examples/Potjans2014
	// Assumptions / limitations:
	// 		As Urgese 2016, we work with a % of the total number of neurons and synapses
	// 		Every neuron in a population receives background activity (8Hz)
	// 		Thalamic current injection is 120Hz and lasts 10 ms as Potjans PyNN example (starts at t_end * 0.7f)
	// 		It does not include the 2x in weight from L4e -> L23e
	// 
	float conn_threshold = 0.0f;		// cap connectivity probability (anything below this value will be set to 0)
	int bg_activity = 8;
	bool thalamic_input = false;
	int thalamic_rate = 120;
	int thalamic_pop = 902;
		
	NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
	cell_params.push_back(c_params);
	// Values from Potjans and Diesmann 2014 paper
	c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
	c_params->v_start_dev = 5.0f;		// std for start voltage (mV) was 5
	c_params->tau_m = 10.0f;			// membrane time constant (ms) 
	c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
	c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
	c_params->v_rest = -65.0f;			// resting potential (mV)
	c_params->v_thres = -50.0f;			// spike threshold (mV) 
	c_params->v_reset = -65.0f;			// voltage after spike (mV) 
	c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
	c_params->i_offset = 0.0f;			// constant current input (nA)
	c_params->c_m = 0.250f;				// cm (nF)
	c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
	c_params->w_exc_dev = 0.0088f;			// nA (PyNN default is nA) was 0.0088
	c_params->w_inh_g = -4.0f;			// g is a multiplier factor with respect to W_EXC
	c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
	c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
	c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms)
	c_params->exc_delay_dev = 0.7f;		// delay on exc synaptic propagation (ms) was 0.7
	c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
	c_params->inh_delay_dev = 0.4f;		// delay on inh synaptic propagation (ms) was 0.4
		
	// populations
	int n_full[8] = {20683,5834,21915, 5479,4850,1065,14395,2948};
	
	Population* l23e = new Population("l2_3e",n_full[0] * n_scale,EXC,c_params); pops.push_back(l23e);
	Population* l23i = new Population("l2_3i",n_full[1] * n_scale,INH,c_params); pops.push_back(l23i);
	Population* l4e = new Population("l4e",n_full[2] * n_scale,EXC,c_params); pops.push_back(l4e);
	Population* l4i = new Population("l4i",n_full[3] * n_scale,INH,c_params); pops.push_back(l4i);
	Population* l5e = new Population("l5e",n_full[4] * n_scale,EXC,c_params); pops.push_back(l5e);
	Population* l5i = new Population("l5i",n_full[5] * n_scale,INH,c_params); pops.push_back(l5i);
	Population* l6e = new Population("l6e",n_full[6] * n_scale,EXC,c_params); pops.push_back(l6e);
	Population* l6i = new Population("l6i",n_full[7] * n_scale,INH,c_params); pops.push_back(l6i);
	
	// connections ([target][source])
	float conn_probs[8][8] = {{0.1009f,  0.1689f, 0.0437f, 0.0818f, 0.0323f, 0.f, 0.0076f, 0.f},
							{0.1346f,   0.1371f, 0.0316f, 0.0515f, 0.0755f, 0.f,     0.0042f, 0.f},
							{0.0077f,   0.0059f, 0.0497f, 0.135f,  0.0067f, 0.0003f, 0.0453f, 0.f},
							{0.0691f,   0.0029f, 0.0794f, 0.1597f, 0.0033f, 0.f,     0.1057f, 0.f},
							{0.1004f,   0.0622f, 0.0505f, 0.0057f, 0.0831f, 0.3726f, 0.0204f, 0.f},
							{0.0548f,   0.0269f, 0.0257f, 0.0022f, 0.06f,   0.3158f, 0.0086f, 0.f},
							{0.0156f,   0.0066f, 0.0211f, 0.0166f, 0.0572f, 0.0197f, 0.0396f, 0.2252f},
							{0.0364f,   0.001f,  0.0034f, 0.0005f, 0.0277f, 0.008f,  0.0658f, 0.1443f}};
		
	// cap connection probabilities
	for(int ii=0; ii < 8; ii++) {
		for(int jj=0; jj < 8; jj++) {
			if(conn_probs[ii][jj] < conn_threshold) conn_probs[ii][jj] = 0;
		}
	}
		
	// in-degree of connectivity
	float k_full[8][8] = {0};
	for(int ii=0; ii < pops.size(); ii++) {
		for(int jj=0; jj < pops.size(); jj++) {
			k_full[ii][jj] = roundf(log(1.f-conn_probs[ii][jj])/log((double)(n_full[ii] * n_full[jj] - 1) / (double)(n_full[ii] * n_full[jj]))) / n_full[ii];
		}
	}
	// build connections
	for(int ii = 0; ii < pops.size(); ii++) {
		for(int jj=0; jj < pops.size(); jj++) {
			Connectivity* c = new RandomConnectivity(pops[ii],pops[jj],k_scale * k_full[jj][ii] * pops[jj]->num_neurons / (float)(pops[ii]->num_neurons * pops[jj]->num_neurons)); //conn_probs[jj][ii] * k_scale);
			conns.push_back(c);
		}
	}
	
	// injections
	int k_ext[8] = {1600,1500,2100,1900,2000,1900,2900,2100};
	for(int ii = 0; ii < pops.size(); ii++) {
		injections.push_back(new PoissonInjector(k_ext[ii] * k_scale,bg_activity,pops[ii],0,t_end,1.0f,dt,k_scale));
	}
	
	// thalamic input 
	if(thalamic_input) {
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[2],t_end * 0.7f,10,0.0983f,dt,k_scale));
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[3],t_end * 0.7f,10,0.0619f,dt,k_scale));
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[6],t_end * 0.7f,10,0.0512f,dt,k_scale));
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[7],t_end * 0.7f,10,0.0196f,dt,k_scale));
	}
	
	// current injection (compensate for k_scale < 1.0)
	if(k_scale < 1.0f) {
		float mean_rate[8] = {0.971f,2.868f,4.746f, 5.396f,8.142f,9.078f,0.991f,7.523f};
		for(int ii=0; ii < pops.size(); ii++) {
			// target pop: pops[ii]
			double x1 = 0;
			// for each source to this target
			for(int jj = 0; jj < pops.size(); jj++) {
				double w;
				if(pops[jj]->type == EXC) w = pops[jj]->cell_params->w_exc_mean;
				else w = pops[jj]->cell_params->w_exc_mean * pops[jj]->cell_params->w_inh_g;
				x1 += w * k_full[ii][jj] * mean_rate[jj]; // weight of connection between pops
			}
			
			x1 += c_params->w_exc_mean * k_ext[ii] * bg_activity; // weight of connection from injector to pop
			double i_ext = 0.5 * (1 - sqrt(k_scale)) * 0.001f * x1; // tau_syn_e * (1-sqrt(k_scale)) * x1
			injections.push_back(new CurrentInjector(1,i_ext,pops[ii],0,t_end,1.0f,dt,1.0f)); 	
		}
	}
	*/
	
	/*
	// USED FOR FRONTIERS IN NEUROINFORMATICS PAPER //
	// CORTICAL MICROCIRCUIT without injectors. To generate activity:
	// v_rest > v_threshold
	// constant input (i_offset)
	float conn_threshold = 0.0f;		// cap connectivity probability (anything below this value will be set to 0)
	
	NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
	cell_params.push_back(c_params);
	// Values from Potjans and Diesmann 2014 paper
	c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
	c_params->v_start_dev = 5.0f;		// std for start voltage (mV) was 5
	c_params->tau_m = 10.0f;			// membrane time constant (ms) 
	c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
	c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
	c_params->v_rest = -45.0f;			// resting potential (mV)
	c_params->v_thres = -50.0f;			// spike threshold (mV) 
	c_params->v_reset = -65.0f;			// voltage after spike (mV) 
	c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
	c_params->i_offset = 0.95f;			// constant current input (nA)
	c_params->c_m = 0.250f;				// cm (nF)
	c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
	c_params->w_exc_dev = 0.0088f;			// nA (PyNN default is nA) was 0.0088
	c_params->w_inh_g = -4.0f;			// g is a multiplier factor with respect to W_EXC
	c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
	c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
	c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms)
	c_params->exc_delay_dev = 0.7f;		// delay on exc synaptic propagation (ms) was 0.7
	c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
	c_params->inh_delay_dev = 0.4f;		// delay on inh synaptic propagation (ms) was 0.4
		
	// populations
	int n_full[8] = {20683,5834,21915, 5479,4850,1065,14395,2948};
	
	Population* l23e = new Population("l2_3e",n_full[0] * n_scale,EXC,c_params); pops.push_back(l23e);
	Population* l23i = new Population("l2_3i",n_full[1] * n_scale,INH,c_params); pops.push_back(l23i);
	Population* l4e = new Population("l4e",n_full[2] * n_scale,EXC,c_params); pops.push_back(l4e);
	Population* l4i = new Population("l4i",n_full[3] * n_scale,INH,c_params); pops.push_back(l4i);
	Population* l5e = new Population("l5e",n_full[4] * n_scale,EXC,c_params); pops.push_back(l5e);
	Population* l5i = new Population("l5i",n_full[5] * n_scale,INH,c_params); pops.push_back(l5i);
	Population* l6e = new Population("l6e",n_full[6] * n_scale,EXC,c_params); pops.push_back(l6e);
	Population* l6i = new Population("l6i",n_full[7] * n_scale,INH,c_params); pops.push_back(l6i);
	
	// connections ([target][source])
	float conn_probs[8][8] = {{0.1009f,  0.1689f, 0.0437f, 0.0818f, 0.0323f, 0.f, 0.0076f, 0.f},
							{0.1346f,   0.1371f, 0.0316f, 0.0515f, 0.0755f, 0.f,     0.0042f, 0.f},
							{0.0077f,   0.0059f, 0.0497f, 0.135f,  0.0067f, 0.0003f, 0.0453f, 0.f},
							{0.0691f,   0.0029f, 0.0794f, 0.1597f, 0.0033f, 0.f,     0.1057f, 0.f},
							{0.1004f,   0.0622f, 0.0505f, 0.0057f, 0.0831f, 0.3726f, 0.0204f, 0.f},
							{0.0548f,   0.0269f, 0.0257f, 0.0022f, 0.06f,   0.3158f, 0.0086f, 0.f},
							{0.0156f,   0.0066f, 0.0211f, 0.0166f, 0.0572f, 0.0197f, 0.0396f, 0.2252f},
							{0.0364f,   0.001f,  0.0034f, 0.0005f, 0.0277f, 0.008f,  0.0658f, 0.1443f}};
	
	// cap connection probabilities
	for(int ii=0; ii < 8; ii++) {
		for(int jj=0; jj < 8; jj++) {
			if(conn_probs[ii][jj] < conn_threshold) conn_probs[ii][jj] = 0;
		}
	}
		
	// in-degree of connectivity
	float k_full[8][8] = {0};
	for(int ii=0; ii < pops.size(); ii++) {
		for(int jj=0; jj < pops.size(); jj++) {
			k_full[ii][jj] = roundf(log(1.f-conn_probs[ii][jj])/(log((double)((n_full[ii] * n_full[jj] - 1)) / (double)((n_full[ii] * n_full[jj]))))) / n_full[ii];
		}
	}
	// build connections
	for(int ii = 0; ii < pops.size(); ii++) {
		for(int jj=0; jj < pops.size(); jj++) {
			Connectivity* c = new RandomConnectivity(pops[ii],pops[jj],k_scale * k_full[jj][ii] * pops[jj]->num_neurons / (float)(pops[ii]->num_neurons * pops[jj]->num_neurons)); //conn_probs[jj][ii] * k_scale);
			conns.push_back(c);
		}
	}
	*/

	/*
	// MULTI AREA MODEL --> 32 areas, each one similar connectivity CM
	// Schmidt 2018 https://github.com/INM-6/multi-area-model
	// Neurons: ~4.13 million * n_scale
	// Synapses: ~24.2 billion * n_scale^2 * k_scale
	// Model only uses type I (intra area) and III (coricocortical) synapses
	// No external input 
		// constant input (i_offset)
	// Cell params from CM
	NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
	cell_params.push_back(c_params);
	// Values from Potjans and Diesmann 2014 paper
	c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
	c_params->v_start_dev = 5.0f;		// std for start voltage (mV) was 5
	c_params->tau_m = 10.0f;			// membrane time constant (ms) 
	c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
	c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
	c_params->v_rest = -65.0f;			// resting potential (mV)
	c_params->v_thres = -50.0f;			// spike threshold (mV) 
	c_params->v_reset = -65.0f;			// voltage after spike (mV) 
	c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
	c_params->i_offset = 0.38f;			// constant current input (nA)
	c_params->c_m = 0.250f;				// cm (nF)
	c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
	c_params->w_exc_dev = 0.0088f;			// nA (PyNN default is nA) was 0.0088
	c_params->w_inh_g = -15.0f;			// g is a multiplier factor with respect to W_EXC was -4
	c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
	c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
	c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms)
	c_params->exc_delay_dev = 0.7f;		// delay on exc synaptic propagation (ms) was 0.7
	c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
	c_params->inh_delay_dev = 0.4f;		// delay on inh synaptic propagation (ms) was 0.4

	// json examples and doc https://en.wikibooks.org/wiki/JsonCpp
	std::ifstream ifs("multiareaModel.json");
	Json::Reader reader;
    Json::Value root;
    reader.parse(ifs, root); // reader can also read strings

	// population sizes in root["neuron_numbers"][AREA][POPULATION]
	// population should be modulated by n_scale
	std::map<std::string,Population*> mappedPops;
	const Json::Value& areas = root["neuron_numbers"]; // array
	std::vector<std::string> areasKeys = areas.getMemberNames();
	int ns = 0;
	for(size_t a=0; a < areasKeys.size(); a++) {
		const std::string& areaKey = areasKeys[a];
		const Json::Value& populations = areas[areaKey];
		std::vector<std::string> populationsKeys = populations.getMemberNames();
		for(size_t p=0; p < populationsKeys.size(); p++) {
			const std::string& popKey = populationsKeys[p];
			float popSize = populations[popKey].asFloat();
			// pop name --> areaKey + popKey
			// pop size --> (int)popSize.asFloat() * n_scale
			// pop type if pop name contains e = EXC, else INH
			// ignore 'total' values
			std::string e_string = "total";
			// discard "total" values and empty population
			if(popKey.find(e_string) != std::string::npos || popSize <= 0) continue;
			e_string = "E";
			std::string pName = areaKey+popKey;
			Population* pop = new Population(pName.c_str(),
										popSize * n_scale,
										popKey.find(e_string) != std::string::npos ? EXC : INH,
										c_params);
			pops.push_back(pop);
			mappedPops[pName] = pop;
			ns += popSize * n_scale;
		}
	}

	// connectivity (number of synapses) in root["synapses"][target_area][target_pop][source_area][source_pop]
	// connectivity should be modulated based on n_scale and k_scale
	const Json::Value& targetAreas = root["synapses"]; // array
	std::vector<std::string> targetAreaKeys = targetAreas.getMemberNames();
	long int syns = 0;
	for(size_t ta=0; ta < targetAreaKeys.size(); ta++) {
		const std::string& taKey = targetAreaKeys[ta];
		const Json::Value& targetPops = targetAreas[taKey];
		if(targetPops.type() != Json::objectValue) continue;
		std::vector<std::string> targetPopsKeys = targetPops.getMemberNames();
		for(size_t tp=0; tp < targetPopsKeys.size(); tp++) {
			const std::string& tpKey = targetPopsKeys[tp];
			const Json::Value& sourceAreas = targetPops[tpKey];
			if(sourceAreas.type() != Json::objectValue) continue;
			std::vector<std::string> sourceAreasKeys = sourceAreas.getMemberNames();
			for(size_t sa=0; sa < sourceAreasKeys.size(); sa++) {
				const std::string& saKey = sourceAreasKeys[sa];
				const Json::Value& sourcePops = sourceAreas[saKey];
				if(sourcePops.type() != Json::objectValue) continue;
				std::vector<std::string> sourcePopsKeys = sourcePops.getMemberNames();
				for(size_t sp=0; sp < sourcePopsKeys.size(); sp++) {
					const std::string& spKey = sourcePopsKeys[sp];
					if(sourcePops[spKey].type() != Json::realValue) continue;
					float number_synapses = sourcePops[spKey].asFloat();
					//printf("[%s][%s][%s][%s] = %f\n",taKey.c_str(),tpKey.c_str(),saKey.c_str(),spKey.c_str(),sourcePops[spKey].asFloat());
					// connectivity: 
						// from taKey + tpKey --> saKey + spKey
						// number synapses --> sourcePops[spKey].asFloat()
					std::string sourcePop = saKey + spKey;
					std::string targetPop = taKey + tpKey;
					if(mappedPops.count(sourcePop) <= 0 || mappedPops.count(targetPop) <= 0) continue; // population was not defined in the population list
					int source_neurons = mappedPops[sourcePop]->num_neurons;
					int target_neurons = mappedPops[targetPop]->num_neurons;
					float ratio = number_synapses / (float)(source_neurons / n_scale * target_neurons / n_scale);
					float p_conn = k_scale * ratio;
					if(p_conn <= 0) {
						continue;
					}
					//printf("%f (%i,%i) --> (%f) final %f\n",number_synapses,source_neurons,target_neurons ,ratio,p_conn);
					// create connection between populations
					Connectivity* c = new RandomConnectivity(
										mappedPops[sourcePop],
										mappedPops[targetPop],
										p_conn
					); 
					conns.push_back(c);
					syns += (long int)((long int)(source_neurons * target_neurons) * p_conn);
				}
			}
		}
	}
	PRINTF("From model --> Number of neurons: %i. Approx synapses: %li\n",ns,syns);
	*/

	
	/*
	// Vogels and Abbot model
	// Parameters based on Brette implementation https://github.com/NeuralEnsemble/PyNN/blob/master/examples/VAbenchmarks.py
	NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
	cell_params.push_back(c_params);
	c_params->v_start_mean = -60.0f;		// mean start voltage (mV)
	c_params->v_start_dev = 5.0f;		// std for start voltage (mV)
	c_params->tau_m = 20.0f;			// membrane time constant (ms) 
	c_params->r_m = 80.0f;				// membrane resistance (MOhm)  is this the right value???
	c_params->tau_ref = 5.0f;			// refractory period after spike (ms) 
	c_params->v_rest = -49.0f;			// resting potential (mV) normal -49
	c_params->v_thres = -50.0f;			// spike threshold (mV) 
	c_params->v_reset = -60.0f;			// voltage after spike (mV) 
	c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
	c_params->i_offset = 0.0f;			// constant current input (nA)
	c_params->c_m = 1000.0f;			// cm (nF)
	c_params->w_exc_mean = 0.0162f;		// nA (PyNN default is nA)
	c_params->w_exc_dev = 0.0f;			// nA (PyNN default is nA)
	c_params->w_inh_g = -5.56f;			// g is a multiplier factor with respect to W_EXC
	c_params->exc_tau_s = 5.0f;			// time constant for exc synapses (ms) 
	c_params->inh_tau_s = 10.0f;		// time constant for inh synapses (ms) 
	c_params->exc_delay_mean = 0.2f;	// delay on exc synaptic propagation (ms)
	c_params->exc_delay_dev = 0.0f;		// delay on exc synaptic propagation (ms)
	c_params->inh_delay_mean = 0.2f;	// delay on inh synaptic propagation (ms)
	c_params->inh_delay_dev = 0.0f;		// delay on inh synaptic propagation (ms)
	
	Population* p_e = new Population("Exc",3200 * n_scale,EXC,c_params);
	Population* p_i = new Population("Inh",800 * n_scale,INH,c_params);
	pops.push_back(p_e);
	pops.push_back(p_i);
	
	Connectivity* e_e = new RandomConnectivity(p_e,p_e,0.02f * k_scale); conns.push_back(e_e);
	Connectivity* e_i = new RandomConnectivity(p_e,p_i,0.02f * k_scale); conns.push_back(e_i);
	Connectivity* i_e = new RandomConnectivity(p_i,p_e,0.02f * k_scale); conns.push_back(i_e);
	Connectivity* i_i = new RandomConnectivity(p_i,p_i,0.02f * k_scale); conns.push_back(i_i);
	
	// for CUBA model, no input (V_REST = -49)
	//PoissonInjector* p_e_inj = new PoissonInjector(20 * k_scale,100,&p_e,0,50,0.01f,dt,k_scale);
	//PoissonInjector* p_i_inj = new PoissonInjector(20 * k_scale,100,&p_i,0,50,0.01f,dt,k_scale);
	//injections.push_back(p_e_inj);
	//injections.push_back(p_i_inj);
	*/
	
	/*
	// CUSTOM MODEL
	float current_injector_p_conn = 0.1f;
	float current_injector_intensity = 0.0f;
	double internal_p_conn = 0.02f;
	
	NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
	cell_params.push_back(c_params);
	c_params->v_start_mean = -60.0f;		// mean start voltage (mV)
	c_params->v_start_dev = 0.0f;		// std for start voltage (mV)
	c_params->tau_m = 20.0f;			// membrane time constant (ms) 
	c_params->r_m = 80.0f;				// membrane resistance (MOhm)  is this the right value???
	c_params->tau_ref = 5.0f;			// refractory period after spike (ms) 
	c_params->v_rest = -49.0f;			// resting potential (mV) normal -49
	c_params->v_thres = -50.0f;			// spike threshold (mV) 
	c_params->v_reset = -60.0f;			// voltage after spike (mV) 
	c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
	c_params->i_offset = 0.0f;			// constant current input (nA)
	c_params->c_m = 1000.0f;			// cm (nF)
	c_params->w_exc_mean = 0.0162f;		// nA (PyNN default is nA)
	c_params->w_exc_dev = 0.0f;			// nA (PyNN default is nA)
	c_params->w_inh_g = -5.56f;			// g is a multiplier factor with respect to W_EXC was -5.56
	c_params->exc_tau_s = 5.0f;			// time constant for exc synapses (ms) 
	c_params->inh_tau_s = 10.0f;		// time constant for inh synapses (ms) 
	c_params->exc_delay_mean = 0.2f;	// delay on exc synaptic propagation (ms)
	c_params->exc_delay_dev = 0.0f;		// delay on exc synaptic propagation (ms)
	c_params->inh_delay_mean = 0.2f;	// delay on inh synaptic propagation (ms)
	c_params->inh_delay_dev = 0.0f;		// delay on inh synaptic propagation (ms)
	
	// each population is a mini VA model
	for(int ii=0; ii < num_processes; ii++) {
		char buf[BUFSIZ];
		snprintf(buf, sizeof(buf), "Exc%d", ii);
		pops.push_back(new Population(buf,200 * n_scale,EXC,c_params));
		snprintf(buf, sizeof(buf), "Inh%d", ii);
		pops.push_back(new Population(buf,50 * n_scale,INH,c_params));
		conns.push_back(new RandomConnectivity(pops[ii*2],pops[ii*2+1],internal_p_conn));
		conns.push_back(new RandomConnectivity(pops[ii*2+1],pops[ii*2],internal_p_conn));
	}
	
	//for(int ii = 0; ii < pops.size(); ii++) {
	//	injections.push_back(new CurrentInjector(1, current_injector_intensity, pops[ii], 0, t_end, current_injector_p_conn,dt,1.0));
	//}
	*/
	
	/*
	// FROM FILE //
	// Cortical Microcircuit (Potjans and Diensmann)
	float conn_threshold = 0.0f;		// cap connectivity probability (anything below this value will be set to 0)
	int bg_activity = 8;
	bool thalamic_input = false;
	int thalamic_rate = 120;
	int thalamic_pop = 902;
	NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
	cell_params.push_back(c_params);
	// Values from Potjans and Diesmann 2014 paper
	c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
	c_params->v_start_dev = 5.0f;		// std for start voltage (mV) was 5
	c_params->tau_m = 10.0f;			// membrane time constant (ms) 
	c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
	c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
	c_params->v_rest = -65.0f;			// resting potential (mV)
	c_params->v_thres = -50.0f;			// spike threshold (mV) 
	c_params->v_reset = -65.0f;			// voltage after spike (mV) 
	c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
	c_params->i_offset = 0.0f;			// constant current input (nA)
	c_params->c_m = 0.250f;				// cm (nF)
	c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
	c_params->w_exc_dev = 0.0088f;			// nA (PyNN default is nA) was 0.0088
	c_params->w_inh_g = -4.0f;			// g is a multiplier factor with respect to W_EXC
	c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
	c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
	c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms)
	c_params->exc_delay_dev = 0.7f;		// delay on exc synaptic propagation (ms) was 0.7
	c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
	c_params->inh_delay_dev = 0.4f;		// delay on inh synaptic propagation (ms) was 0.4
	
	// ONLY SUPPORTS MODELS LOADED FROM SINGLE FILES (i.e. only one conns object in the list)
	load_model_from_python_file("resources/model.txt",c_params,&pops,&conns);
	float conn_probs[8][8] = {{0.1009f,  0.1689f, 0.0437f, 0.0818f, 0.0323f, 0.f, 0.0076f, 0.f},
							{0.1346f,   0.1371f, 0.0316f, 0.0515f, 0.0755f, 0.f,     0.0042f, 0.f},
							{0.0077f,   0.0059f, 0.0497f, 0.135f,  0.0067f, 0.0003f, 0.0453f, 0.f},
							{0.0691f,   0.0029f, 0.0794f, 0.1597f, 0.0033f, 0.f,     0.1057f, 0.f},
							{0.1004f,   0.0622f, 0.0505f, 0.0057f, 0.0831f, 0.3726f, 0.0204f, 0.f},
							{0.0548f,   0.0269f, 0.0257f, 0.0022f, 0.06f,   0.3158f, 0.0086f, 0.f},
							{0.0156f,   0.0066f, 0.0211f, 0.0166f, 0.0572f, 0.0197f, 0.0396f, 0.2252f},
							{0.0364f,   0.001f,  0.0034f, 0.0005f, 0.0277f, 0.008f,  0.0658f, 0.1443f}};
	// populations
	int n_full[8] = {20683,5834,21915, 5479,4850,1065,14395,2948};
	// in-degree of connectivity
	float k_full[8][8] = {0};
	for(int ii=0; ii < pops.size(); ii++) {
		for(int jj=0; jj < pops.size(); jj++) {
			k_full[ii][jj] = roundf(log(1.f-conn_probs[ii][jj])/log((double)(n_full[ii] * n_full[jj] - 1) / (double)(n_full[ii] * n_full[jj]))) / n_full[ii];
		}
	}
	// injections
	int k_ext[8] = {1600,1500,2100,1900,2000,1900,2900,2100};
	for(int ii = 0; ii < pops.size(); ii++) {
		injections.push_back(new PoissonInjector(k_ext[ii] * k_scale,bg_activity,pops[ii],0,t_end,1.0f,dt,k_scale));
	}
	// thalamic input 
	if(thalamic_input) {
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[2],t_end * 0.7f,10,0.0983f,dt,k_scale));
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[3],t_end * 0.7f,10,0.0619f,dt,k_scale));
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[6],t_end * 0.7f,10,0.0512f,dt,k_scale));
		injections.push_back(new PoissonInjector(thalamic_pop * k_scale,thalamic_rate,pops[7],t_end * 0.7f,10,0.0196f,dt,k_scale));
	}
	
	// current injection (compensate for k_scale < 1.0)
	if(k_scale < 1.0f) {
		float mean_rate[8] = {0.971f,2.868f,4.746f, 5.396f,8.142f,9.078f,0.991f,7.523f};
		for(int ii=0; ii < pops.size(); ii++) {
			// target pop: pops[ii]
			double x1 = 0;
			// for each source to this target
			for(int jj = 0; jj < pops.size(); jj++) {
				double w;
				if(pops[jj]->type == EXC) w = pops[jj]->cell_params->w_exc_mean;
				else w = pops[jj]->cell_params->w_exc_mean * pops[jj]->cell_params->w_inh_g;
				x1 += w * k_full[ii][jj] * mean_rate[jj]; // weight of connection between pops
			}
			
			x1 += c_params->w_exc_mean * k_ext[ii] * bg_activity; // weight of connection from injector to pop
			double i_ext = 0.5 * (1 - sqrt(k_scale)) * 0.001f * x1; // tau_syn_e * (1-sqrt(k_scale)) * x1
			injections.push_back(new CurrentInjector(1,i_ext,pops[ii],0,t_end,1.0f,dt,1.0f)); 	
		}
	}
	*/


	/*
	// FROM FILE//
	// Fruit Fly brain //
	NeuronParams* c_params = (NeuronParams*)malloc(sizeof(NeuronParams));
	cell_params.push_back(c_params);
	c_params->v_start_mean = -58.0f;		// mean start voltage (mV)
	c_params->v_start_dev = 0.0f;		// std for start voltage (mV) was 5
	c_params->tau_m = 10.0f;			// membrane time constant (ms) 
	c_params->r_m = 40.0f;				// membrane resistance (MOhm) 
	c_params->tau_ref = 2.0f;			// refractory period after spike (ms) 
	c_params->v_rest = -49.0f;			// resting potential (mV) was -65
	c_params->v_thres = -50.0f;			// spike threshold (mV) 
	c_params->v_reset = -65.0f;			// voltage after spike (mV) 
	c_params->v_spike = 20.0f;			// nominal spiking potential (mV) (only for drawing purposes) 
	c_params->i_offset = 0.0f;			// constant current input (nA)
	c_params->c_m = 0.250f;				// cm (nF)
	c_params->w_exc_mean = 0.0878f;		// nA (PyNN default is nA)
	c_params->w_exc_dev = 0.0f;			// nA (PyNN default is nA) was 0.0088
	c_params->w_inh_g = -4.0f;			// g is a multiplier factor with respect to W_EXC
	c_params->exc_tau_s = 0.5f;			// time constant for exc synapses (ms) 
	c_params->inh_tau_s = 0.5f;		// time constant for inh synapses (ms) 
	c_params->exc_delay_mean = 1.5f;	// delay on exc synaptic propagation (ms) was 1.5
	c_params->exc_delay_dev = 0.0f;		// delay on exc synaptic propagation (ms) was 0.7
	c_params->inh_delay_mean = 0.8f;	// delay on inh synaptic propagation (ms)
	c_params->inh_delay_dev = 0.0f;		// delay on inh synaptic propagation (ms) was 0.4
	
	// ONLY SUPPORTS MODELS LOADED FROM SINGLE FILES (i.e. only one conns object in the list)
	load_model_from_python_file("resources/v1_2_reduced.txt",c_params,&pops,&conns);
	*/
	
	
	// END OF USER INPUT //
	
	// ASSERT SIMULATION PARAMETERS //
	// propagation time step will be reset to the min synaptic delay and multiple of dt
	float propagation_step = dt;
	if(cell_params.size() == 0 || cell_params[0]->exc_delay_dev > 0 || cell_params[0]->inh_delay_dev > 0) {
		// TODO: find the lowest random syn delay; ensure randomness produces multiples of dt
		PRINTF("%i: Random synaptic delays detected. Setting propagation time step to 1 (potentially expensive)\n",process_id);
	} else {
		// synaptic delay must be a multiple of dt
		float exc_rem = remainderf(cell_params[0]->exc_delay_mean,dt);
		float inh_rem = remainderf(cell_params[0]->inh_delay_mean,dt);
		if(fabs(exc_rem) > 0 || fabs(inh_rem) > 0) {
			PRINTF("%i: Automatically adjusting synaptic delays\n",process_id);
			PRINTF("%i: EXC: %f || INH: %f\n",process_id,cell_params[0]->exc_delay_mean,cell_params[0]->inh_delay_mean);
			cell_params[0]->exc_delay_mean -= exc_rem;
			cell_params[0]->inh_delay_mean -= inh_rem;
		}  
		propagation_step = min(cell_params[0]->exc_delay_mean,cell_params[0]->inh_delay_mean); 
	}
	// ensure propagation time step is within bounds (2^5)
	// time-driven propagation step must be a multiple of dt (ensured since syn delays are multiples of dt)
	int propagation_time_step = min(propagation_step / dt + 0.5f,pow(2,5)); // (in number of sim timesteps)
	
	PRINTF("%i: Propagation step: %i\n",process_id,propagation_time_step);
	// Max neuron ID (2^27-1) (if using RMA communicator, cannot use the last int value
	int total_neurons = 0;
	for(int ii=0; ii < pops.size(); ii++) {
			total_neurons += pops[ii]->num_neurons;
	}
	if(total_neurons >= pow(2,27)) { // limit imposed by multiplexing propagation messages (5 bits spike time, 27 bits neuron id)
		printf("%i: Simulator only supports 2^27 distinct neurons (current: %i)\n",process_id,total_neurons);
		MPI_Finalize();
		return 0;
	}	
	// END ASSERTION //

	// Create network (connections and injections)
	// network
	int population_size = 0;
	for(int ii=0; ii < pops.size(); ii++) {
		population_size += pops[ii]->num_neurons;
	}

	if(!MODEL_LOADED_FROM_FILE) {
		// if the model is loaded from file, do not randomise IDs
		// if model is created here, randomise IDs
		// this is to avoid issues with IDs particularly when detecting which neuron belongs to which population
		// generate randomised neuron ids
		// It is not a problem when creating networks here as the actual connection pattern is not created
		// until later; but if loaded from file, there is a strong link between pre and post synaptic neuron ids
		std::vector<int> used_ids(population_size,-1);
		for(int ii=0; ii < population_size; ii++) {
			used_ids[ii] = ii;
		}
		std::random_shuffle(used_ids.begin(),used_ids.end());
		for(int ii=0; ii < pops.size(); ii++) {
			pops[ii]->neuron_ids.clear();
			pops[ii]->neuron_ids.swap(pops[ii]->neuron_ids);
			pops[ii]->neuron_ids.resize(pops[ii]->num_neurons);
			for(int jj=0; jj < pops[ii]->num_neurons; jj++) {
				int number = used_ids.back();
				pops[ii]->neuron_ids[jj] = number;
				used_ids.pop_back();
			}
		}
	}


	// injections
	srand(random_seed); // ensure base population with same seed has same injectors connectivity
	for(int ii=0; ii < injections.size(); ii++) {
		injections[ii]->create_connections(process_id);
	}
	
	// Create partition object to hold connections and partitioning information across processes
	// Assign neurons to partitions
	/*Partitioning* partition;
	if(strcmp(partitioning,"graphPartitioning") == 0) {
		PRINTF("%i: Partitioning: graph partition (METIS)\n",process_id);
		partition = new GraphPartitionPartitioning(&pops,population_size);
	} else if(strcmp(partitioning,"randomBalanced") == 0) {  
		PRINTF("%i: Partitioning: random-balanced\n",process_id);
		partition = new RandomBalancedPartitioning(&pops,population_size);
	} else if(strcmp(partitioning,"parallelPartitioning") == 0) {  
		PRINTF("%i: Partitioning: parallel partitioning (PARMETIS)\n",process_id);
		partition = new ParallelPartitionPartitioning(&pops,population_size);
	} else if(strcmp(partitioning,"hypergraphPartitioning") == 0) {  
		PRINTF("%i: Partitioning: hypergraph partitioning (Zoltan)\n",process_id);
		partition = new HypergraphPartitioning(&pops,population_size);
	} else if(strcmp(partitioning,"hierPartitioning") == 0) {  
		PRINTF("%i: Partitioning: hypergraph hierarchical partitioning (Zoltan)\n",process_id);
		partition = new HIERPartitioning(&pops,population_size);
	} else if(strcmp(partitioning,"custom") == 0) {  
		PRINTF("%i: Partitioning: custom\n",process_id);
		if(custom_partitioning == NULL) partition = new RandomPartitioning(&pops,population_size);
		else partition = new CustomPartitioning(&pops,population_size,custom_partitioning);
	} else if(strcmp(partitioning,"praw") == 0) {  
		PRINTF("%i: Partitioning: PRAW\n",process_id);
		partition = new PRAWFilePartitioning(&pops,population_size,comm_bandwidth_matrix_file);
	} else if(strcmp(partitioning,"zoltanFile") == 0) {  
		PRINTF("%i: Partitioning: Zoltan from file\n",process_id);
		partition = new ZoltanFilePartitioning(&pops,population_size);
	} else if(strcmp(partitioning,"roundrobin") == 0) {  
		PRINTF("%i: Partitioning: Round robin\n",process_id);
		partition = new RoundRobinPartitioning(&pops,population_size);
	} else {
		PRINTF("%i: Partitioning: random\n",process_id);
		partition = new RandomPartitioning(&pops,population_size);
	}*/
	
	
	// load neuron activity from file (if provided)
	std::vector<int> previous_neuron_activity(population_size,1);
	if(use_weight_file) {
		PRINTF("%i: Using previous neuron activity in partitioning...\n",process_id);
		if(process_id == MASTER_NODE) {
			// load from file
			std::ifstream inputFile(w_file);
			if (inputFile) {        
				int value;
				// read the elements in the file into a vector  
				int counter=0;
				while ( inputFile >> value ) {
					previous_neuron_activity[counter] = value;
					counter++;
				}
				inputFile.close();
			} else {
				printf("%i: Neuron activity file %s not found!\n",process_id,argv[5]);
			}
			// send values to other nodes
			for(int ii=0; ii < num_processes; ii++) {
				if(ii == process_id) continue;
				MPI_Send(&(previous_neuron_activity[0]),previous_neuron_activity.size(), MPI_INT,ii,0,MPI_COMM_WORLD);
			}
		} else {
			int* buf = (int*)malloc(sizeof(int) * population_size);
			MPI_Recv(buf, population_size, MPI_INT, MASTER_NODE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int ii=0; ii < population_size; ii++) {
				previous_neuron_activity[ii] = buf[ii];
			}
			free(buf);
		}
		
	}

	// construct Model from connectivity and population objects
	// to save memory, construct only in one process per node, 
	// then distribute complete model to neighbours
	int exc = 0, inh = 0;
	model.population_size = population_size;
	model.interconnections_size = (int*)calloc(model.population_size,sizeof(int));
	model.interconnections = (int**)malloc(sizeof(int*) * model.population_size);
	if(process_id % node_size == 0) {
		// get the number of connections per pre synaptic neuron
		for(int ii=0; ii < conns.size(); ii++) {
			PRINTF("%i: Creating connection set %i out of %lu\n",process_id,ii+1,conns.size());
			if(!MODEL_LOADED_FROM_FILE) conns[ii]->generate_connectivity(model.population_size);
			for(int jj=0; jj < conns[ii]->connections.size(); jj++) {
				int con_length = conns[ii]->connections[jj].size();
				int nid = conns[ii]->from->neuron_ids[jj];
				model.interconnections_size[nid] += con_length;
			}
			
		}
		// initialise pre-synaptic neuron data structures
		int* lastId = (int*)calloc(model.population_size,sizeof(int));
		for(int ii=0; ii < model.population_size; ii++) {
			model.interconnections[ii] = (int*)malloc(sizeof(int) * model.interconnections_size[ii]);
		}
		// assign post-synaptic neuron ids (and weights) to connection lists
		for(int ii=0; ii < conns.size(); ii++) {
			for(int jj=0; jj < conns[ii]->connections.size(); jj++) {
				int con_length = conns[ii]->connections[jj].size();
				int nid = conns[ii]->from->neuron_ids[jj];
				for(int cc=0; cc < con_length; cc++) {
					if(conns[ii]->from->type == EXC) {
						model.interconnections[nid][lastId[nid]] = conns[ii]->connections[jj][cc];
						exc++;
					} else {
						model.interconnections[nid][lastId[nid]] = -conns[ii]->connections[jj][cc];
						inh++;
					}
					lastId[nid] += 1;
				}
			}
			
			// clean up connectivity object
			delete conns[ii];
		}
		conns.clear();
		free(lastId);

		// send model to other processes
		for(int ii=0; ii < node_size; ii++) {
			int target_process = ii + process_id / node_size * node_size;
			if(target_process == process_id) continue;
			// exc and inh connection numbers
			MPI_Send(&exc,1,MPI_INT,target_process,0,MPI_COMM_WORLD);
			MPI_Send(&inh,1,MPI_INT,target_process,0,MPI_COMM_WORLD);
			// model.interconnections_size
			MPI_Send(model.interconnections_size,model.population_size,MPI_INT,target_process,0,MPI_COMM_WORLD);
			for(int jj=0; jj < model.population_size; jj++) {
				MPI_Send(model.interconnections[jj],model.interconnections_size[jj],MPI_INT,target_process,jj,MPI_COMM_WORLD);
			}
		}
		
	} else {
		// receive model from central process in node
		int central_process = process_id / node_size * node_size;
		// exc and inh connections
		MPI_Recv(&exc,1,MPI_INT,central_process,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&inh,1,MPI_INT,central_process,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		// model.interconnections_size
		MPI_Recv(model.interconnections_size,model.population_size,MPI_INT,central_process,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		// model.interconnections
		for(int ii=0; ii < model.population_size; ii++) {
			model.interconnections[ii] = (int*)malloc(sizeof(int) * model.interconnections_size[ii]);
		}
		for(int ii=0; ii < model.population_size; ii++) {
			MPI_Recv(model.interconnections[ii],model.interconnections_size[ii],MPI_INT,central_process,ii,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		// clean up conns structure
		for(int ii=0; ii < conns.size(); ii++) {
			delete conns[ii];
		}
		conns.clear();
	}
	
	
	/*int exc = 0, inh = 0;
	model.population_size = population_size;
	// get the number of connections per pre synaptic neuron
	model.interconnections_size = (int*)calloc(model.population_size,sizeof(int));
	for(int ii=0; ii < conns.size(); ii++) {
		PRINTF("%i: Creating connection set %i out of %lu\n",process_id,ii+1,conns.size());
		if(!MODEL_LOADED_FROM_FILE) conns[ii]->generate_connectivity(model.population_size);
		for(int jj=0; jj < conns[ii]->connections.size(); jj++) {
			int con_length = conns[ii]->connections[jj].size();
			int nid = conns[ii]->from->neuron_ids[jj];
			model.interconnections_size[nid] += con_length;
		}
	}
	// initialise pre-synaptic neuron data structures
	int* lastId = (int*)calloc(model.population_size,sizeof(int));
	model.interconnections = (int**)malloc(sizeof(int*) * model.population_size);
	for(int ii=0; ii < model.population_size; ii++) {
		model.interconnections[ii] = (int*)malloc(sizeof(int) * model.interconnections_size[ii]);
	}
	// assign post-synaptic neuron ids (and weights) to connection lists
	for(int ii=0; ii < conns.size(); ii++) {
		for(int jj=0; jj < conns[ii]->connections.size(); jj++) {
			int con_length = conns[ii]->connections[jj].size();
			int nid = conns[ii]->from->neuron_ids[jj];
			for(int cc=0; cc < con_length; cc++) {
				if(conns[ii]->from->type == EXC) {
					model.interconnections[nid][lastId[nid]] = conns[ii]->connections[jj][cc];
					exc++;
				} else {
					model.interconnections[nid][lastId[nid]] = -conns[ii]->connections[jj][cc];
					inh++;
				}
				lastId[nid] += 1;
			}
		}
		
		// clean up connectivity object
		delete conns[ii];
	}
	conns.clear();
	free(lastId);
	*/


	model.total_connections = exc + inh;
	PRINTF("%i: Neurons: %lu, Total connections: %i (%i exc, %i inh)\n", process_id,model.population_size, model.total_connections,exc,inh);
	
	if(model.store_in_file && process_id == MASTER_NODE) {
		// store model in file
		PRINTF("%i: Storing model in file resources/model.txt...\n",process_id);
		FILE *fp = fopen("resources/model.txt", "w+");
		
		// write header
		fprintf(fp,"%lu",model.population_size);
		for(int ii=0; ii < pops.size(); ii++) {
			if(pops[ii]->type == EXC)
				fprintf(fp," %i",pops[ii]->num_neurons);
			else
				fprintf(fp," -%i",pops[ii]->num_neurons);
		}
		fprintf(fp,"\n");
		// write connectivity per neuron id
		for(int ii=0; ii < model.population_size; ii++) {
			// find what pop ii belongs to
			for(int p=0; p < pops.size(); p++) {
				if(pops[p]->is_in_population(ii)) {
					fprintf(fp,"%i",p);
					break;
				}
			}
			// write all connections from ii
			for(int jj=0; jj < model.interconnections_size[ii]; jj++) {
				fprintf(fp," %i",abs(model.interconnections[ii][jj]));
			}
			fprintf(fp,"\n");
		}
		
		fclose(fp);
	
	}

	
	// run as many simulations as permutations of partitioning methods x comm patterns
	// all share the same model (built previously)
	int iterations = partitioning_methods.size() * comm_pattern_methods.size();
	for(int coms=0; coms < comm_pattern_methods.size(); coms++) {
		////////////////////////////
		// START OF PERMUTATION LOOP
		////////////////////////////
		for(int pats=0; pats < partitioning_methods.size();pats++) {
			double global_timer = MPI_Wtime();
			double timer;
			
			int iteration = coms * comm_pattern_methods.size() + pats;
			const char* comm_method = comm_pattern_methods[coms].c_str();
			const char* part_method = partitioning_methods[pats].c_str();

			srand(random_seed);
			randgauss(0,0,true);

			// Create partition object to hold connections and partitioning information across processes
			// Assign neurons to partitions
			Partitioning* partition;
			if(strcmp(part_method,"graphPartitioning") == 0) {
				PRINTF("%i: Partitioning: graph partition (METIS)\n",process_id);
				partition = new GraphPartitionPartitioning(&pops,population_size);
			} else if(strcmp(part_method,"randomBalanced") == 0) {  
				PRINTF("%i: Partitioning: random-balanced\n",process_id);
				partition = new RandomBalancedPartitioning(&pops,population_size);
			} else if(strcmp(part_method,"parallelPartitioning") == 0) {  
				PRINTF("%i: Partitioning: parallel partitioning (PARMETIS)\n",process_id);
				partition = new ParallelPartitionPartitioning(&pops,population_size);
			} else if(strcmp(part_method,"hypergraphPartitioning") == 0) {  
				PRINTF("%i: Partitioning: hypergraph partitioning (Zoltan)\n",process_id);
				partition = new HypergraphPartitioning(&pops,population_size);
			} else if(strcmp(part_method,"hierPartitioning") == 0) {  
				PRINTF("%i: Partitioning: hypergraph hierarchical partitioning (Zoltan)\n",process_id);
				partition = new HIERPartitioning(&pops,population_size);
			} else if(strcmp(part_method,"custom") == 0) {  
				PRINTF("%i: Partitioning: custom\n",process_id);
				if(custom_partitioning == NULL) partition = new RandomPartitioning(&pops,population_size);
				else partition = new CustomPartitioning(&pops,population_size,custom_partitioning);
			} else if(strcmp(part_method,"prawS_without") == 0) {  
				PRINTF("%i: Partitioning: sequential PRAW\n",process_id);
				partition = new PRAWFilePartitioning(&pops,population_size,comm_bandwidth_matrix_file,false,false);
			} else if(strcmp(part_method,"prawP_without") == 0) {  
				PRINTF("%i: Partitioning: parallel PRAW\n",process_id);
				partition = new PRAWFilePartitioning(&pops,population_size,comm_bandwidth_matrix_file,true,false);
			} else if(strcmp(part_method,"prawS") == 0) {  
				PRINTF("%i: Partitioning: sequential PRAW\n",process_id);
				partition = new PRAWFilePartitioning(&pops,population_size,comm_bandwidth_matrix_file,false,true);
			} else if(strcmp(part_method,"prawP") == 0) {  
				PRINTF("%i: Partitioning: parallel PRAW\n",process_id);
				partition = new PRAWFilePartitioning(&pops,population_size,comm_bandwidth_matrix_file,true,false);
			} else if(strcmp(part_method,"zoltanFile") == 0) {  
				PRINTF("%i: Partitioning: Zoltan from file\n",process_id);
				partition = new ZoltanFilePartitioning(&pops,population_size,comm_bandwidth_matrix_file);
			} else if(strcmp(part_method,"roundrobin") == 0) {  
				PRINTF("%i: Partitioning: Round robin\n",process_id);
				partition = new RoundRobinPartitioning(&pops,population_size);
			} else {
				PRINTF("%i: Partitioning: random\n",process_id);
				partition = new RandomPartitioning(&pops,population_size);
			}
			
			partition->perform_partitioning(&model,num_processes,process_id,&previous_neuron_activity);

			// fast access info
			idx_t* parts = partition->partitioning;
			int pop_size = model.population_size;
			
			// Communication strategy
			Communicator* communicator;	
			if(strcmp(comm_method,"subscriber") == 0) {
				PRINTF("%i: Comm pattern: subscriber\n",process_id);
				communicator = new SubscriberCommunicator(process_id, num_processes, MASTER_NODE,model.population_size);
			} else if (strcmp(comm_method,"direct") == 0) {
				communicator = new DirectCommunicator(process_id, num_processes, MASTER_NODE, model.population_size);
				PRINTF("%i: Comm pattern: DirectP2P\n",process_id);
			} else if(strcmp(comm_method,"rma") == 0) {
				communicator = new RMACommunicator(process_id, num_processes, MASTER_NODE,model.population_size);
				PRINTF("%i: Comm pattern: RMA\n",process_id);
			} else if(strcmp(comm_method,"nbx") == 0) {
				communicator = new NBXCommunicator(process_id, num_processes, MASTER_NODE,model.population_size);
				PRINTF("%i: Comm pattern: NBX\n",process_id);
			} else if(strcmp(comm_method,"pcx") == 0) {
				communicator = new PCXCommunicator(process_id, num_processes, MASTER_NODE,model.population_size);
				PRINTF("%i: Comm pattern: PCX\n",process_id);
			} else if(strcmp(comm_method,"rsx") == 0) {
				communicator = new RSXCommunicator(process_id, num_processes, MASTER_NODE,model.population_size);
				PRINTF("%i: Comm pattern: RSX\n",process_id);
			} else if(strcmp(comm_method,"pex") == 0) {
				communicator = new PEXCommunicator(process_id, num_processes, MASTER_NODE,model.population_size);
				PRINTF("%i: Comm pattern: PEX\n",process_id);
			} else {
				communicator = new AllGatherCommunicator(process_id, num_processes, MASTER_NODE,model.population_size);
				PRINTF("%i: Comm pattern: All gather\n",process_id);
			}
			
			
			// Cost of partition cut (total)
			// total computational load imbalance (as well as neuron and synaptic balances)
			// calculated as max imbalance / average imbalance
			// takes into account 1 unit per hosted neuron and 1 unit per incoming synapse to hosted neurons
			// store the graph connectivity (list of connected partitions)
			int* total_local_synapses = (int*)calloc(num_processes,sizeof(int));
			int total_cut = 0;
			int* neurons_in_partition = (int*)calloc(num_processes,sizeof(int));
			std::vector<std::set<int> > partition_connectivity(num_processes);
			for(int ii = 0; ii < model.population_size; ii++) {
				neurons_in_partition[parts[ii]] += 1;
				for(int jj=0; jj < model.interconnections_size[ii]; jj++) {
					int part_number = parts[abs(model.interconnections[ii][jj])];
					total_local_synapses[part_number] += 1;
					if(parts[ii] != part_number) {
						total_cut++;
						partition_connectivity[parts[ii]].insert(parts[abs(model.interconnections[ii][jj])]);
					}
				}	
			}
			double synapse_workload = 0;
			int max_synapse_workload = 0;
			double total_workload = 0;
			int max_total_workload = 0;
			double neuron_workload = 0;
			int max_neuron_workload = 0;
			for(int ii=0; ii < num_processes; ii++) {
				total_workload += neurons_in_partition[ii] + total_local_synapses[ii];
				if(neurons_in_partition[ii] + total_local_synapses[ii] > max_total_workload)
					max_total_workload = neurons_in_partition[ii] + total_local_synapses[ii];
				neuron_workload += neurons_in_partition[ii];
				if(neurons_in_partition[ii] > max_neuron_workload) 
					max_neuron_workload = neurons_in_partition[ii];
				synapse_workload += total_local_synapses[ii];
				if(max_synapse_workload < total_local_synapses[ii])
					max_synapse_workload = total_local_synapses[ii];
			}
			neuron_workload /= num_processes;
			total_workload /= num_processes;
			synapse_workload /= num_processes;
			double neuron_comp_imbalance = neuron_workload > 0 ? max_neuron_workload / neuron_workload: 0;
			double total_comp_imbalance = total_workload > 0 ? max_total_workload / total_workload : 0;
			double synapse_comp_imbalance = synapse_workload > 0 ? max_synapse_workload / synapse_workload : 0;

			PRINTF("%i: Synapses cut: %i\n",process_id,total_cut);
			PRINTF("%i: Neuronal computational imbalance: %f\n",process_id,neuron_comp_imbalance);
			PRINTF("%i: Synaptic computational imbalance: %f\n",process_id,synapse_comp_imbalance);
			PRINTF("%i: Total workload imbalance %f (%f)\n",process_id,total_comp_imbalance,total_workload);

			free(total_local_synapses);
			
		#ifdef PARTITION_CONNECTIVITY_GRAPH
			// store partition connectivity graph in a file
			if(process_id == MASTER_NODE) { 
				std::string filename = sim_name;
				char str_int[16];
				filename += "_";
				filename += part_method;
				filename += "_";
				filename += comm_method;
				filename += "_partition_graph_";
				sprintf(str_int,"%i",num_processes);
				filename +=  str_int;

				FILE *fp = fopen(filename.c_str(), "w+");
				fprintf(fp,"%i\n",num_processes);
				for(int ii=0; ii < num_processes; ii++) {
					fprintf(fp,"%lu ",partition_connectivity[ii].size());
					for (std::set<int>::iterator it = partition_connectivity[ii].begin(); it != partition_connectivity[ii].end(); ++it) {
						fprintf(fp,"%i ",*it);
					}
					fprintf(fp,"\n");
				}
				fclose(fp);
			}
		#endif
			
			// Potentially UNSAFE: since same seed, all processes do the partition equally
			
			//////
			// create data structures (Neurons, Synapses and local CurrentInjectors) per process
			//////
			
			//synchronise random generator
			srand(random_seed);
			
			// initialise neuron population and current injector
			PRINTF("%i: Initialising neurons and connections...\n",process_id);
			// globalNeuronRefs --> array of global IDs to quickly locate local neurons on each partition
			Neuron** globalNeuronRefs = NULL;
			// list of local neurons only (starts full then will be trimmed
			std::vector<Neuron> localNeurons;
			int local_neurons = 0;
			if(!model.null_compute) {
				globalNeuronRefs = (Neuron**)calloc(model.population_size,sizeof(Neuron*));
				localNeurons.resize(neurons_in_partition[process_id]);
				for(int ii=0; ii < pops.size(); ii++) {
					for(int jj=0; jj < pops[ii]->neuron_ids.size(); jj++) {
						int neuronId = pops[ii]->neuron_ids[jj];
						// we do the init voltage for all so the entire network shares initial voltage values 
						// irrespective of the number of partitions
						double init_v = randgauss(pops[ii]->cell_params->v_start_mean,pops[ii]->cell_params->v_start_dev);	// random start potential
						init_v -= remainderf(init_v,0.1f); // round off voltage to one decimal (helps with float precision errors)
						if(parts[neuronId] == process_id) {
							// neuron is local
							localNeurons[local_neurons].init(neuronId,pops[ii]->cell_params,init_v);
							globalNeuronRefs[neuronId] = &localNeurons[local_neurons];
							local_neurons++;
						}
					}
				}
			}
			free(neurons_in_partition);
			
			PRINTF("%i: Local neurons: %i\n",process_id,local_neurons);
			// UNSAFE: random() may not be equal across processes since each process called random() (in Neuron constructor) an unequal amount of times
			
			
			// create injection reservoires
			PRINTF("%i: Culling injection reservoires\n",process_id);
			for(int ii=0; ii < injections.size(); ii++) {
				injections[ii]->cull_to_local_post_syn_neurons(partition->partitioning,process_id);
			}
			
			// create target synapses for neuron (only if target is local process)
			PRINTF("%i: Creating synapse counts\n",process_id);
			int local_synapses = 0;
			std::vector<Synapse> localSynapses;
			std::vector<int> globalSynapse_xadj;
			std::vector<int> globalSynapse_adjncy;
			// count number of postsynaptic targets for each neuron
			//std::vector<int> num_synapses(model.population_size,0);
			for(int ii=0; ii < model.population_size; ii++) {
				for(int jj=0; jj < model.interconnections_size[ii]; jj++) {
					// ii -> presynaptic neuron; interconnections[ii][jj] --> postsynaptic neuron
					// notify communicator
					if(previous_neuron_activity[ii] > 0) // use previous activity info to determine if a link is going to become active
						communicator->from_to_connection(ii,abs(model.interconnections[ii][jj]),parts);
					if(parts[abs(model.interconnections[ii][jj])] == process_id) { 
						//num_synapses[ii] += 1;
						local_synapses++;
					}
				}
			}
			
			if(!model.null_compute) {	
				PRINTF("%i: Creating synapse objects and references\n",process_id);
				localSynapses.resize(local_synapses); // all local synapses
				int syn_count;
				if(local_synapses > 0) { // if there are no local synapses, skip assigning memory for them
					syn_count = 0;
					bool check_cell_params = cell_params.size() > 1;
					globalSynapse_xadj.resize(pop_size+1);
					for(int ii=0; ii < pop_size; ii++) {
						int conn_size = model.interconnections_size[ii];
						int jj=0;
						globalSynapse_xadj[ii] = syn_count;
						//for(int* con = model.interconnections[ii], jj = 0; jj < conn_size; jj++,con++) {
						for(int p=0; p < conn_size; p++) {
							int postneuron = abs(model.interconnections[ii][p]);
							// ii -> presynaptic neuron; interconnections[ii][jj] --> postsynaptic neuron
							// which population interconnections[ii][jj] belongs to
							NeuronParams* c = NULL;
							if(check_cell_params) {
								for(int kk=0; kk < pops.size(); kk++) {
									if(pops[kk]->is_in_population(postneuron)) {
										c = pops[kk]->cell_params;
										break;
									}
								}
							} else c = cell_params[0];
							// weight and delay of the synapse calculated irrespective of whether synapse will reside in local process
							// ensures identical network init across seeded sims
							float weight = randgauss(c->w_exc_mean, c->w_exc_dev);
							weight = weight < 0 ? 0 : weight;
							float delay = 0;
							bool isExcitatory = model.interconnections[ii][p] > 0;
							if(!isExcitatory) {
								weight *= c->w_inh_g;
								delay = randgauss(c->inh_delay_mean, c->inh_delay_dev);
							} else {
								delay = randgauss(c->exc_delay_mean, c->exc_delay_dev);
							}

							delay = delay < 0 ? 0 : delay;

							if(parts[postneuron] == process_id) { 
								localSynapses[syn_count].init(ii, globalNeuronRefs[postneuron],isExcitatory,dt,c,k_scale,weight,delay);
								syn_count++;
							}
							
						}
					}
					globalSynapse_xadj[pop_size] = syn_count;
					
					
					std::vector<int> indices(local_synapses);
					for(int ii=0; ii < local_synapses; ii++) {
						indices[ii] = ii;
					}
					std::sort(indices.begin(),indices.end(),IdxCompare(localSynapses));
					// sorting by the postsynaptic neuron id is important
					// to make the synapse update step faster (utilising mem cache)
					std::sort(localSynapses.begin(),localSynapses.end()); // (fly brain) sort causes sim with seed 1516885998 to have different spike count (234138) as expected (234144)
					
					globalSynapse_adjncy.resize(local_synapses);
					for(int ii = 0; ii < local_synapses; ii++) {
						globalSynapse_adjncy[indices[ii]] = ii;
					}
				}

				PRINTF("%i: Local synapses: %i\n",process_id,local_synapses);
			}
			
			// clear unused model connection data
			// relevant process-process dependencies are stored in communicator->neuron_connection_to
			PRINTF("%i: Clearing local data\n",process_id);
			if(iteration >= iterations-1){
				for(int ii=0; ii < model.population_size; ii++) {
					free(model.interconnections[ii]);		
				}
				free(model.interconnections_size);
				free(model.interconnections);
			}
			//simulate
			PRINTF("%i: Preparing start simulation...\n",process_id);
			long int numSpikes = 0;
			long int remote_spikes = 0;
			int simSize = (int)(t_end / dt);

			
			//propagation objects
			std::vector<unsigned int> propagate;
			std::vector<int> fired;
			std::vector<unsigned int> received;
			
			//monitoring objects
			double* computation_times = (double*)calloc(simSize,sizeof(double));
			if(computation_times == NULL) exit_error();
			double* propagate_times = (double*)calloc(simSize,sizeof(double));
			if(propagate_times == NULL) exit_error();
			double* sync_time = (double*)calloc(simSize / propagation_time_step,sizeof(double));
			if(sync_time == NULL) exit_error();
			double* idle_time = (double*)calloc(simSize / propagation_time_step,sizeof(double));
			if(idle_time == NULL) exit_error();
			std::vector<std::vector<int> > neuron_activity; // monitor local neuron spike activity
		#if defined(RECORD_STATISTICS) || defined(RECORD_ACTIVITY)
			neuron_activity.resize(local_neurons);
			// reserve memory to avoid reallocation during sim
			for(int ii=0; ii < local_neurons; ii++) {
				neuron_activity[ii].reserve(simSize);
			}
		#endif
			
			// trace where processes spend time on each timestep
			int num_tracing_events = 7;
			int num_traces = simSize * num_tracing_events + 1;
			double* trace_timings = NULL;
		#ifdef RECORD_PROCESS_TRACE
			trace_timings = (double*)calloc(num_traces,sizeof(double));
			if(trace_timings == NULL) exit_error();
		#endif
			
			
		#ifdef API_PROFILING
			SCOREP_RECORDING_ON();
		#endif


			// final communicator setup
			if(communicator != NULL) communicator->final_setup();
			
			PRINTF("%i: Waiting for other nodes\n",process_id);
			
			// synchronise all processes
			MPI_Barrier(MPI_COMM_WORLD);
			
			// synchronise mpi time (with respect to MASTER_NODE time)
			double process_time_correction = -MPI_Wtime();
			PRINTF("%i: Start simulation. %f time correction\n",process_id,process_time_correction);
			
			double init_time;
			timer = MPI_Wtime();
			init_time = timer - global_timer;
			global_timer = timer;
			
			// loop variables
			int ii, jj, kk, t;
			float current_t;
			
			for(t=0; t < simSize; t++) {
			
				timer = MPI_Wtime();
		#ifdef RECORD_PROCESS_TRACE
				trace_timings[t * num_tracing_events] = timer + process_time_correction;
		#endif

				current_t = dt * t;
				
		#ifdef API_PROFILING
				SCOREP_USER_REGION_DEFINE(updateInjector_region);
				SCOREP_USER_REGION_BEGIN(updateInjector_region,"updateInjectors", SCOREP_USER_REGION_TYPE_COMMON);
		#endif
				if(!model.null_compute) {
					// update current from injectors
					for(ii=0; ii < injections.size(); ii++) {
						injections[ii]->update(current_t);
						for(jj=0; jj < injections[ii]->size; jj++) { // SLOW!
							for(kk=0; kk < injections[ii]->post_syn_neuron_idx[jj].size(); kk++) {
								globalNeuronRefs[injections[ii]->post_syn_neuron_idx[jj][kk]]->i_e += (injections[ii]->i_e[jj]);
							}
						}
					}
				}
				
		#ifdef API_PROFILING
				SCOREP_USER_REGION_END(updateInjector_region);
		#endif

		#ifdef RECORD_PROCESS_TRACE
				trace_timings[t * num_tracing_events+1] = MPI_Wtime() + process_time_correction;
		#endif

		#ifdef API_PROFILING
				SCOREP_USER_REGION_DEFINE(updateLocalSynapses_region);
				SCOREP_USER_REGION_BEGIN(updateLocalSynapses_region,"updateLocalSynapses", SCOREP_USER_REGION_TYPE_COMMON);
		#endif
				
				// update local synapses 
				// storing local neuron references in a contiguous array saves each process
				// from having to comb the entire g_synapses array
				if(!model.null_compute) {
					Synapse* syn;
					for(syn = &localSynapses[0], ii=0; ii < local_synapses; ii++) {
						syn->update_current(current_t);
						syn++;
					}
				}
				
				
		#ifdef API_PROFILING
				SCOREP_USER_REGION_END(updateLocalSynapses_region);
		#endif

		#ifdef RECORD_PROCESS_TRACE
				trace_timings[t * num_tracing_events+2] = MPI_Wtime() + process_time_correction;
		#endif

		#ifdef API_PROFILING
				SCOREP_USER_REGION_DEFINE(updateLocalNeuron_region);
				SCOREP_USER_REGION_BEGIN(updateLocalNeuron_region,"updateLocalNeuron", SCOREP_USER_REGION_TYPE_COMMON);
		#endif
				
				// solve neuron state
				fired.clear();
				if(!model.null_compute) {
					// storing local neuron references in a contiguous array saves each process
					// from having to comb the entire neurons array and testing if local (most of which won't be local)
					Neuron* neuron;
					for(neuron = &localNeurons[0], ii=0; ii < local_neurons; ii++) {
						if(neuron->update_potential(dt,current_t)) {
							numSpikes++;
							fired.push_back(neuron->id);
		#if defined(RECORD_STATISTICS) || defined(RECORD_ACTIVITY)
							neuron_activity[ii].push_back(t);
		#endif
						}
						neuron++;
					}
				} else {
					// no compute simulation, add all local neurons to propagation list
					for(ii=0; ii < model.population_size; ii++) {
						if(parts[ii] == process_id) {
							fired.push_back(ii);
							numSpikes++;
						}
					}
				}
				
		#ifdef API_PROFILING
				SCOREP_USER_REGION_END(updateLocalNeuron_region);
		#endif

				computation_times[t] = MPI_Wtime() - timer;
				timer = MPI_Wtime();
		#ifdef RECORD_PROCESS_TRACE
				trace_timings[t * num_tracing_events+3] = timer + process_time_correction;
		#endif

		#ifdef API_PROFILING
				SCOREP_USER_REGION_DEFINE(propagate_region);
				SCOREP_USER_REGION_BEGIN(propagate_region,"propagate", SCOREP_USER_REGION_TYPE_COMMON);
		#endif	

				// gather spikes
				received.clear();
				for(ii=0; ii < fired.size(); ii++) {
					int preneuron = fired[ii];
					if(!model.null_compute) {
						//first go through local targets
						if(globalSynapse_xadj[preneuron] < globalSynapse_xadj[preneuron+1]) {
							received.push_back(encodeSpike(0,preneuron));
						}
					}
					// check if pre-synaptic firing neuron has any listener partitions
					if(communicator->neuron_connection_to[preneuron].size() > 0) {
						remote_spikes++;
						// add to propagate the time of spiking and neuron id
						// how many timesteps until the next propagation_time_step
						propagate.push_back(encodeSpike(propagation_time_step - (t % propagation_time_step) - 1,preneuron));
					}
				}

				// only propagate every propagation_time_step
				if((t+1) % propagation_time_step == 0) {
					// even though communicator->send_receive already coordinates processes
					// processes that may not communicate on one time step may need to do so in the next
					// if no barrier, then they may go out of sync
					double sub_t;
		#ifdef RECORD_PROCESS_TRACE
					trace_timings[t * num_tracing_events+4] = MPI_Wtime() + process_time_correction;
		#endif
		#ifdef MEASURE_IDLE_TIME
					// Calling MPI_Barrier can improve performance by reducing network contention
					sub_t = MPI_Wtime();
					MPI_Barrier(MPI_COMM_WORLD);
					idle_time[(t+1) / propagation_time_step - 1] = MPI_Wtime() - sub_t;
					
		#endif
					// disseminate spikes across processes
					sub_t = MPI_Wtime();
		#ifdef RECORD_PROCESS_TRACE
					trace_timings[t * num_tracing_events+5] = sub_t + process_time_correction;
		#endif
					//communicator->send_receive(&propagate, parts, &(model.connections), &received);
					communicator->send_receive(&propagate, &received);
					sync_time[(t+1) / propagation_time_step - 1] = MPI_Wtime() - sub_t;
					propagate.clear();

		#ifdef MEASURE_IDLE_TIME
					// this barrier is necessary to avoid imbalance from send_receive to spill over to idle time
					// Or not. MPI_Barrier syncs order of events, but does not ensure process are in timed synced
					// Processes reaching barrier first send message to all others and wait
					// Last processes to reach would acknowledge and continue; and send message to all others
					// The processes arriving to the barrier first will wait for the messages to arrive from slower processes
					// By the time they receive the messages, the system is likely to be imbalanced again
					// Suggested: Use Simple Time Protocol to sync in time
					// https://www.techopedia.com/definition/4539/simple-network-time-protocol-sntp
					//MPI_Barrier(MPI_COMM_WORLD);
		#endif
					
				} else {
		#ifdef RECORD_PROCESS_TRACE
					// still record tracing timings
					double timing = MPI_Wtime();
					trace_timings[t * num_tracing_events+4] = timing + process_time_correction;
					trace_timings[t * num_tracing_events+5] = timing + process_time_correction;
		#endif
				}
				
					
		#ifdef API_PROFILING
				SCOREP_USER_REGION_END(propagate_region);
		#endif
				
				propagate_times[t] = MPI_Wtime() - timer;// measure time the process was propagating spikes
				
				timer = MPI_Wtime();
		#ifdef RECORD_PROCESS_TRACE
				trace_timings[t * num_tracing_events+6] = timer + process_time_correction;
		#endif

		#ifdef API_PROFILING
				SCOREP_USER_REGION_DEFINE(updateSynapse_region);
				SCOREP_USER_REGION_BEGIN(updateSynapse_region,"updateSynapse", SCOREP_USER_REGION_TYPE_COMMON);
		#endif
				if(!model.null_compute) {
					// update synapses
					for(ii=0; ii < received.size(); ii++) {
						// decode neuron id and spike timing
						int preneuron;
						int spike_timing;
						decodeSpike(received[ii],&spike_timing,&preneuron);
						for(jj=globalSynapse_xadj[preneuron]; jj < globalSynapse_xadj[preneuron+1]; jj++) {
							localSynapses[globalSynapse_adjncy[jj]].spike(spike_timing);
						}
					}
				}
				
				
				
		#ifdef API_PROFILING
				SCOREP_USER_REGION_END(updateSynapse_region);
		#endif	
				computation_times[t] += MPI_Wtime() - timer; // add update synapse computing times to neuron and synapse update 		
				
				
			}
		#ifdef RECORD_PROCESS_TRACE
			trace_timings[simSize * num_tracing_events] = MPI_Wtime() + process_time_correction;
		#endif

			double sim_time;
			timer = MPI_Wtime();
			sim_time = timer - global_timer;
			global_timer = timer;
			
		#ifdef API_PROFILING
			SCOREP_USER_REGION_DEFINE(gather_region);
			SCOREP_USER_REGION_BEGIN(gather_region,"gatherInfo", SCOREP_USER_REGION_TYPE_COMMON);
		#endif

			// gather info from processes (spikes and timings)
			PRINTF("End simulation. Start gathering results...\n");

			long int total_spikes;
			MPI_Allreduce(&numSpikes, &total_spikes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
			long int total_remote_spikes;
			MPI_Allreduce(&remote_spikes, &total_remote_spikes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
			
			double comptime = 0;
			double proptime = 0;
			double idletime = 0;
			double synctime = 0;
			for(ii = 0; ii < simSize; ii++) {
				comptime += computation_times[ii];
				proptime += propagate_times[ii];
				if(ii < simSize / propagation_time_step) {
					idletime += idle_time[ii];
					synctime += sync_time[ii];
				}
			}
			
			std::vector<std::vector<double> > total_traces;
		#ifdef RECORD_PROCESS_TRACE
			// share trace timings per process //
			if(process_id != MASTER_NODE) {
				// send traces ti MASTER_NODE
				MPI_Send(trace_timings,num_traces, MPI_DOUBLE,MASTER_NODE,process_id,MPI_COMM_WORLD);
			} else {
				total_traces.resize(num_processes);
				// receive one trace per process
				double* buf = (double*)malloc(sizeof(double) * num_traces);
				for(ii=0; ii < num_processes; ii++) {
					if(ii == MASTER_NODE) {
						// collect MASTER_NODE data too
						for(jj=0; jj < num_traces; jj++) {
							total_traces[ii].push_back(trace_timings[jj]);
						}
					} else {
						MPI_Recv(buf, num_traces, MPI_DOUBLE, ii, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						for(jj=0; jj < num_traces; jj++) {
							total_traces[ii].push_back(buf[jj]);
						}
					}
				}
				free(buf);
			}
			////
		#endif

			/* variability in timings throughout simulation */
			double comp_variance = 0;
			double idle_variance = 0;
			double sync_variance = 0;
			for(ii = 0; ii < simSize; ii++) {
				comp_variance += (computation_times[ii] - comptime/simSize) * (computation_times[ii] - comptime/simSize);
				if(ii < simSize / propagation_time_step) {
					idle_variance += (idle_time[ii] - (idletime/(simSize/propagation_time_step))) * (idle_time[ii] - (idletime/(simSize/propagation_time_step)));
					sync_variance += (sync_time[ii] - (synctime/(simSize/propagation_time_step))) * (sync_time[ii] - (synctime/(simSize/propagation_time_step)));
				}
			}
			comp_variance = sqrt(comp_variance / simSize);
			idle_variance = sqrt(idle_variance / (simSize / propagation_time_step));
			sync_variance = sqrt(sync_variance / (simSize / propagation_time_step));

			// Runtime neighbours during simulation
			long int* runtime_neighbours = (long int*)malloc(sizeof(long int) * num_processes);
			MPI_Gather(&communicator->runtime_neighbours,1,MPI_LONG,runtime_neighbours,1,MPI_LONG,MASTER_NODE,MPI_COMM_WORLD);
			int* non_empty_comm_steps = (int*)malloc(sizeof(int) * num_processes);
			MPI_Gather(&communicator->non_empty_comm_step,1,MPI_INT,non_empty_comm_steps,1,MPI_INT,MASTER_NODE,MPI_COMM_WORLD);
			long int total_runtime_neighbours = 0;
			int total_non_empty_comm_steps = 0;
			for(ii=0; ii < num_processes; ii++){
				total_runtime_neighbours += runtime_neighbours[ii];
				total_non_empty_comm_steps += non_empty_comm_steps[ii];
			}
			float average_runtime_neighbours = (float)total_runtime_neighbours / num_processes / (num_processes-1) / communicator->comm_step;
			float average_non_empty_runtime_neighbours = total_non_empty_comm_steps == 0 ? 0 : (float)total_runtime_neighbours / (num_processes-1) / total_non_empty_comm_steps;
			
			int total_sync_group_size;
			int sync_group_size = communicator->get_sync_group_size();
			MPI_Allreduce(&sync_group_size, &total_sync_group_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			double avg_sync_group_size = (double)total_sync_group_size / num_processes / (num_processes-1);
			double* sim_times = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&sim_time,1,MPI_DOUBLE,sim_times,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			double* comp_times = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&comptime,1,MPI_DOUBLE,comp_times,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			double* idle_times = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&idletime,1,MPI_DOUBLE,idle_times,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			//int* total_local_synapses = (int*)malloc(sizeof(int)*num_processes);
			//MPI_Gather(&local_synapses,1,MPI_INT,total_local_synapses,1,MPI_INT,MASTER_NODE,MPI_COMM_WORLD);
			//double total_synapse_comp_imbalance;
			//MPI_Allreduce(&synapse_comp_imbalance, &total_synapse_comp_imbalance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			double* p_times = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&proptime,1,MPI_DOUBLE,p_times,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			long int* comm_sent = (long int*)malloc(sizeof(long int) * num_processes);
			MPI_Gather(&communicator->communication_sent,1,MPI_LONG,comm_sent,1,MPI_LONG,MASTER_NODE,MPI_COMM_WORLD);
			double* s_times = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&synctime,1,MPI_DOUBLE,s_times,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			long int messages_sent;
			MPI_Allreduce(&communicator->num_messages, &messages_sent, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
			long int total_interprocess_spikes;
			MPI_Allreduce(&(communicator->spikes_sent), &total_interprocess_spikes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
			double* comp_variances = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&comp_variance,1,MPI_DOUBLE,comp_variances,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			double* idle_variances = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&idle_variance,1,MPI_DOUBLE,idle_variances,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			double* sync_variances = (double*)malloc(sizeof(double) * num_processes);
			MPI_Gather(&sync_variance,1,MPI_DOUBLE,sync_variances,1,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
			
			std::vector<std::vector<int> > global_activity;
			
		#if defined(RECORD_STATISTICS) || defined(RECORD_ACTIVITY)
			// Share neuron activity record with master node
			if(process_id != MASTER_NODE) {
				// send neuron_activity to MASTER_NODE
				for(ii = 0; ii < neuron_activity.size(); ii++) {
					if(neuron_activity[ii].size() == 0) continue;
					int gid = localNeurons[ii].id;
					MPI_Send(&neuron_activity[ii][0],neuron_activity[ii].size(),MPI_INT,MASTER_NODE,gid,MPI_COMM_WORLD);
				}
				// finished message
				MPI_Send(MPI_BOTTOM,0,MPI_INT,MASTER_NODE,model.population_size,MPI_COMM_WORLD);
			} else {
				global_activity.resize(model.population_size);
				// add local activity
				for(ii=0; ii < local_neurons; ii++) {
					int gid = localNeurons[ii].id;
					for(jj=0; jj < neuron_activity[ii].size(); jj++) {
						global_activity[gid].push_back(neuron_activity[ii][jj]);
					}
				}
				// receive neuron_activity from all nodes
				for(ii=0; ii < num_processes; ii++) {
					if(ii == MASTER_NODE) continue;
					bool listen = true;
					MPI_Status status;
					int count;
					do {
						MPI_Probe(ii,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
						MPI_Get_count(&status, MPI_INT, &count);
						if(status.MPI_TAG == model.population_size) {
							// termination message, receive and ignore
							int ignore;
							MPI_Recv(&ignore, 0, MPI_INT, ii, status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
							listen = false;
						} else {
							// receive activity series for neuron status.MPI_TAG
							int* buffer = (int*)malloc(sizeof(int) * count);
							MPI_Recv(buffer, count, MPI_INT, ii, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							for(int ss=0; ss < count; ss++) {
								global_activity[status.MPI_TAG].push_back(buffer[ss]);
							}
							free(buffer);
						}
						
					} while(listen);
				}
			}
			
		#endif	

		#if defined(ADVANCED_COMM_STATS) || defined(ADVANCED_COMM_STATS_MATRIX_ONLY)
			// store advanced communication stats (weight of each message)
			// EXPERIMENT_NAME_comm_sizes_NUM_PROCESSES will contain sizes of all messages (from all processes)
			// EXPERIMENT_NAME_comm_matrix_NUM_PROCESSES will contain table with total comm between any two processes
			if(process_id == MASTER_NODE) {
				std::string filename = sim_name;
				filename += "_";
				filename += part_method;
				filename += "_";
				filename += comm_method;
				char str_int[16];
				sprintf(str_int,"%i",num_processes);
				filename += "_comm_sizes_";
				filename +=  str_int;
				FILE *fp;
		#ifndef ADVANCED_COMM_STATS_MATRIX_ONLY
				fp = fopen(filename.c_str(), "w+");
		#endif
				MPI_Status status;
				int count;
				// store values for MASTER_NODE first
				int* current_comm_matrix = (int*)calloc(num_processes,sizeof(int));
				for(ii=0; ii < communicator->messages_sent.size(); ii++) {
					int total = 0;
					for(int p=0; p < communicator->messages_sent[ii].size(); p++) {
						total += communicator->messages_sent[ii][p];
		#ifndef ADVANCED_COMM_STATS_MATRIX_ONLY
						fprintf(fp,"%i ",communicator->messages_sent[ii][p]);
		#endif
					}
					current_comm_matrix[ii] = total;
					communicator->messages_sent[ii].clear();
					communicator->messages_sent[ii].swap(communicator->messages_sent[ii]);
				}
		#ifndef ADVANCED_COMM_STATS_MATRIX_ONLY
				// collect message values from other nodes
				for(ii=0; ii < num_processes; ii++) {
					if(ii == MASTER_NODE) continue;
					MPI_Probe(ii,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Get_count(&status, MPI_INT, &count);
					int* buffer = (int*)malloc(sizeof(int) * count);
					MPI_Recv(buffer, count, MPI_INT, ii, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(int ss=0; ss < count; ss++) {
						// write each value to a file
						fprintf(fp,"%i ",buffer[ss]);
					}
					free(buffer);
				}
				fprintf(fp,"\n");
				fclose(fp);

				MPI_Barrier(MPI_COMM_WORLD); // make sure first communication batch is completed
		#endif
				// collect total values from other processes (comm matrix)
				filename = sim_name;
				filename += "_";
				filename += part_method;
				filename += "_";
				filename += comm_method;
				filename += "_comm_matrix_";
				filename +=  str_int;
				fp = fopen(filename.c_str(), "w+");
				// store MASTER_NODE local values first
				for(ii=0; ii < num_processes; ii++) {
					fprintf(fp,"%i ",current_comm_matrix[ii]);
				}
				fprintf(fp,"\n");
				// collect message values from other nodes
				for(ii=0; ii < num_processes; ii++) {
					if(ii == MASTER_NODE) continue;
					MPI_Probe(ii,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Get_count(&status, MPI_INT, &count);
					int* buffer = (int*)malloc(sizeof(int) * count);
					MPI_Recv(buffer, count, MPI_INT, ii, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(int ss=0; ss < count; ss++) {
						// write each value to a file
						fprintf(fp,"%i ",buffer[ss]);
					}
					fprintf(fp,"\n");
					free(buffer);
				}

				free(current_comm_matrix);
				fclose(fp);
			} else {
				std::vector<int> flatten;
				int* current_comm_matrix = (int*)calloc(num_processes,sizeof(int));
				for(ii=0; ii < communicator->messages_sent.size(); ii++) {
					int total = 0;
					for(int p=0; p < communicator->messages_sent[ii].size(); p++) {
						total += communicator->messages_sent[ii][p];
						flatten.push_back(communicator->messages_sent[ii][p]);
					}
					current_comm_matrix[ii] = total;
					communicator->messages_sent[ii].clear();
					communicator->messages_sent[ii].swap(communicator->messages_sent[ii]);
				}
		#ifndef ADVANCED_COMM_STATS_MATRIX_ONLY
				int* buf;
				if(flatten.size() > 0)
					buf = &flatten[0];
				MPI_Send(buf,flatten.size(),MPI_INT,MASTER_NODE,process_id,MPI_COMM_WORLD);

				MPI_Barrier(MPI_COMM_WORLD); // make sure first communication batch is completed
		#endif
				MPI_Send(current_comm_matrix,num_processes,MPI_INT,MASTER_NODE,process_id,MPI_COMM_WORLD);
			}
		#endif
			
			double gather_time;
			if(process_id == MASTER_NODE) {
				timer = MPI_Wtime();
				gather_time = timer - global_timer;
				global_timer = timer;
			}
		#ifdef API_PROFILING
			SCOREP_USER_REGION_END(gather_region);
			SCOREP_RECORDING_OFF();
		#endif

			// report
			if(process_id == MASTER_NODE) {
				std::string filename = sim_name;
				filename += "_";
				filename += part_method;
				filename += "_";
				filename += comm_method;
				char str_int[16];
				
		#if defined(RECORD_PROCESSES_ACTIVITY) || defined(RECORD_PROCESS_TRACE)
				// record process average timings + process traces
				// record timings for computation, idle and synchronisation
				// one file per process: EXPERIMENTNAME_NPROCESSORS_PROCESSNUMBER
				sprintf(str_int,"%i",num_processes);
				filename += "_";
				filename +=  str_int;
				for(ii=0; ii < num_processes; ii++) {
					std::string fprocessor = filename;
					fprocessor += "_";
					sprintf(str_int,"%i",ii);
					fprocessor += str_int;
					bool fileexists = access(fprocessor.c_str(), F_OK) != -1;
					FILE *fp = fopen(fprocessor.c_str(), "ab+");
					if(fp == NULL) {
						printf("Error when storing processor results into file\n");
					} else {
						if(!fileexists) // file does not exist, add header
							fprintf(fp,"%s,%s,%s,%s,%s,%s,%s,%s\n","Comp time","Idle time","Sync time","Comp deviation","Idle deviation","Sync deviation","Runtime neighbours","Non-empty runtime neighbours");
						float non_empty_cs = non_empty_comm_steps[ii] == 0 ? 0 : (float)runtime_neighbours[ii]/(num_processes-1)/non_empty_comm_steps[ii];
						fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f\n",comp_times[ii],idle_times[ii],s_times[ii],comp_variances[ii],idle_variances[ii],sync_variances[ii],(float)runtime_neighbours[ii]/(num_processes-1)/communicator->comm_step,non_empty_cs);
					}
					fclose(fp);
		#ifdef RECORD_PROCESS_TRACE
					std::string fprocess_traces = fprocessor;
					fprocess_traces += "_trace";
					fp = fopen(fprocess_traces.c_str(), "ab+");
					if(fp == NULL) {
						printf("Error when storing processor results into file\n");
					} else {
						for(jj=0; jj < total_traces[ii].size()-1; jj++) {
							fprintf(fp,"%.9f,",total_traces[ii][jj]);
						}
						fprintf(fp,"%.9f\n",total_traces[ii][total_traces.size()-1]);
					}
					fclose(fp);
		#endif
				}		
				
		#endif		

				double spikeRate = (1000.0f * total_spikes) / t_end;
				printf("Spikes: %li, (%li neurons), \nAverage firing rate per neuron: %f Hz\n",total_spikes,model.population_size,spikeRate/model.population_size);
				printf("Init time: %fs; Gather info time: %fs\n",init_time,gather_time);
				printf("Computation times:\n");
				double total_comp = 0;
				//double min_prop = RAND_MAX;
				//double max_prop = 0;
				double total_prop = 0;
				double total_sync_time = 0;
				double max_sim_time = 0;
				double total_idle_time = 0;
				
				
				for(ii = 0; ii < num_processes; ii++){
					PRINTF("%i: %fs (%f%%); %fs propagate; %fs sync\n",ii,comp_times[ii],comp_times[ii]/sim_time*100,p_times[ii],s_times[ii]);
					//if(p_times[ii] > max_prop)
					//	max_prop = p_times[ii];
					//if(p_times[ii] < min_prop)
					//	min_prop = p_times[ii];
					total_prop += p_times[ii];
					total_sync_time += s_times[ii];
					total_comp += comp_times[ii];
					if(max_sim_time < sim_times[ii])
						max_sim_time = sim_times[ii];
					total_idle_time += idle_times[ii];
				}
				// calculate variances
				double prop_variance = 0;
				double comp_variance = 0;
				double sync_variance = 0;
				double idle_variance = 0;
				double avg_comp = total_comp/num_processes;
				double avg_prop = total_prop/num_processes;
				double avg_sync = total_sync_time/num_processes;
				double avg_idle = total_idle_time/num_processes;
				for(ii = 0; ii < num_processes; ii++) {
					comp_variance += (comp_times[ii]-avg_comp) * (comp_times[ii] - avg_comp);
					prop_variance += (p_times[ii]-avg_prop) * (p_times[ii] - avg_prop);
					sync_variance += (s_times[ii]-avg_sync) * (s_times[ii] - avg_sync);
					idle_variance += (idle_times[ii]-avg_idle) * (idle_times[ii] - avg_idle);
				}
				comp_variance = sqrt(comp_variance / num_processes);
				prop_variance = sqrt(prop_variance / num_processes);
				sync_variance = sqrt(sync_variance / num_processes);
				idle_variance = sqrt(idle_variance / num_processes);

				printf("Deviation in computation: %f\n",comp_variance);
				printf("Deviation in propagation: %f\n",prop_variance);
				printf("Deviation in synchronisation (data exchange): %f\n",sync_variance);
				printf("Deviation in idle (implicit sync): %f\n",idle_variance);
				printf("Theoretical - messages sent (bytes):\n");
				long int total_sent = 0;
				for(ii = 0; ii < num_processes; ii++){
					total_sent += comm_sent[ii];
					PRINTF("%i: %lu bytes\n",ii,comm_sent[ii]);
				}
				printf("Total sim time: %f\n",max_sim_time);
				printf("Total communication: %lu messages, %lu bytes\n",messages_sent,total_sent);
				printf("Total sync time (per process): %f\n",total_sync_time/num_processes);
				printf("Total idle time (per process): %f\n",total_idle_time/num_processes);
				printf("Total computational time (per process): %f\n",avg_comp);
				// record global sim results
				filename = sim_name;
				filename += "_";
				filename += part_method;
				filename += "_";
				filename += comm_method;
				sprintf(str_int,"%i",num_processes);
				filename += "__";
				filename +=  str_int;
				bool fileexists = access(filename.c_str(), F_OK) != -1;
				FILE *fp = fopen(filename.c_str(), "ab+");
				double tcm = (double)total_interprocess_spikes/(double)total_spikes;
				if(fp == NULL) {
					printf("Error when storing results into file\n");
				} else {
					/*if(!fileexists) // file does not exist, add header
						fprintf(fp,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","Build time","Sim time (per process)","Comp time (per process)","Comp variance","Propagation time (per process)","Prop variance","Sync time (per process)","Idle time (per process)","Bytes sent (minimum),Messages sent","Random seed","Edgecut","Total connections","Edgecut ratio","Neuronal computation imbalance","Synaptic computation imbalance","Spikes sent","Total spikes","Remote spikes","Total comm cost","Comm cost ratio","Avg sync group size","Average runtime neighbours","Average non-empty runtime neighbours");
					fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%lu,%lu,%i,%i,%i,%f,%f,%f,%lu,%i,%i,%f,%f,%f,%f,%f\n",init_time,max_sim_time,avg_comp,comp_variance,avg_prop,prop_variance,total_sync_time/num_processes,total_idle_time/num_processes,total_sent,messages_sent,random_seed,total_cut,model.total_connections,(float)total_cut/model.total_connections,neuron_comp_imbalance,total_synapse_comp_imbalance,total_interprocess_spikes,total_spikes,total_remote_spikes,tcm,tcm / (num_processes-1),avg_sync_group_size,average_runtime_neighbours,average_non_empty_runtime_neighbours);
					*/
					if(!fileexists) // file does not exist, add header
						fprintf(fp,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","Build time","Sim time (per process)","Comp time (per process)","Comp deviation","Sync deviation","Idle deviation","Sync time (per process)","Idle time (per process)","Bytes sent (minimum),Messages sent","Random seed","Edgecut","Total connections","Edgecut ratio","Total computation imbalance","Neuronal computation imbalance","Synaptic computation imbalance","Spikes sent","Total spikes","Remote spikes","Total comm cost","Comm cost ratio","Avg sync group size","Average runtime neighbours","Average non-empty runtime neighbours");
					fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%lu,%lu,%i,%i,%i,%f,%f,%f,%f,%lu,%li,%li,%f,%f,%f,%f,%f\n",init_time,max_sim_time,avg_comp,comp_variance,sync_variance,idle_variance,total_sync_time/num_processes,total_idle_time/num_processes,total_sent,messages_sent,random_seed,total_cut,model.total_connections,(float)total_cut/model.total_connections,total_comp_imbalance,neuron_comp_imbalance,synapse_comp_imbalance,total_interprocess_spikes,total_spikes,total_remote_spikes,tcm,tcm / (num_processes-1),avg_sync_group_size,average_runtime_neighbours,average_non_empty_runtime_neighbours);

				}
				fclose(fp);
		#ifdef RECORD_ACTIVITY		
				// neuron_activity contains the spike activity (timings) of all neurons
				// store neuron_activity in file
				std::string neuron_activity_filename = filename;
				neuron_activity_filename += "_neuron_activity";
				FILE *fneuron_activity = fopen(neuron_activity_filename.c_str(), "wb"); // all spikes
				for(int ii=0; ii < global_activity.size(); ii++) {
					fprintf(fneuron_activity,"%lu ",global_activity[ii].size()); // ISI
				}
		#endif
		#ifdef RECORD_STATISTICS
				
				// store statistics
				std::string intervals_filename = filename;
				intervals_filename += "_ISI";
				std::string spikes_filename = filename;
				spikes_filename += "_spikes";
				std::string isicv_filename = filename;
				isicv_filename += "_ISICV";
				FILE *fintervals = fopen(intervals_filename.c_str(), "wb"); // all neurons ISI
				FILE *fspikes = fopen(spikes_filename.c_str(), "wb"); // all neurons spikes
				FILE** fisicv_pop = (FILE**)malloc(sizeof(FILE*) * pops.size()); // ISICV per population
				FILE *fisicv = fopen(isicv_filename.c_str(), "wb"); // all neurons ISI CV
				FILE** fintervals_pop = (FILE**)malloc(sizeof(FILE*) * pops.size()); // neurons in population ISI
				FILE** fspikes_pop = (FILE**)malloc(sizeof(FILE*) * pops.size()); // neurons in population spikes
				for(int ii=0; ii < pops.size(); ii++) {
					std::string pop_intervals_filename = intervals_filename;
					pop_intervals_filename += "_";
					pop_intervals_filename += pops[ii]->name;
					fintervals_pop[ii] = fopen(pop_intervals_filename.c_str(), "wb");
					std::string pop_spikes_filename = spikes_filename;
					pop_spikes_filename += "_";
					pop_spikes_filename += pops[ii]->name;
					fspikes_pop[ii] = fopen(pop_spikes_filename.c_str(), "wb");
					std::string pop_isicv_filename = isicv_filename;
					pop_isicv_filename += "_";
					pop_isicv_filename += pops[ii]->name;
					fisicv_pop[ii] = fopen(pop_isicv_filename.c_str(), "wb");
				}
				// check all FILE streams are not NULL!
				// store neuron spikes
				std::vector<std::vector<float> > neuron_isi(model.population_size,std::vector<float>(0));
				std::vector<int> spikes_pop(pops.size(),0);
				for(int ii=0; ii < global_activity.size(); ii++) {
					int spikes = 0;
					float lastSpike = 0;
					int pop_belonging = 0;
					for(int pp=0; pp < pops.size(); pp++) {
						if(pops[pp]->is_in_population(ii)) { 
						//if(ii >= pops[pp]->first_idx && ii < pops[pp]->first_idx + pops[pp]->num_neurons) {
							pop_belonging = pp;
							break;
						}
					}
					for(jj=0; jj < global_activity[ii].size(); jj++) {
						spikes++;
						spikes_pop[pop_belonging] += 1;
						if(lastSpike > 0) {
							neuron_isi[ii].push_back((global_activity[ii][jj] - lastSpike) * dt);
							fprintf(fintervals,"%f ",(global_activity[ii][jj] - lastSpike) * dt); // ISI
							// per pop
							fprintf(fintervals_pop[pop_belonging],"%f ",(global_activity[ii][jj] - lastSpike) * dt);
						}
						lastSpike = global_activity[ii][jj];
						
					}
					fprintf(fspikes,"%i ",spikes); // activity (spikes) of that neuron
					
					// record TOTAL ISI coefficient of variation for that cell (std / mean)
					if(neuron_isi[ii].size() > 1) {
						double total = std::accumulate(neuron_isi[ii].begin(),neuron_isi[ii].end(),0.0);
						double mean = total / neuron_isi[ii].size();
						std::vector<double> diff(neuron_isi[ii].size());
						std::transform(neuron_isi[ii].begin(), neuron_isi[ii].end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
						double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
						double stdev = sqrt(sq_sum / neuron_isi[ii].size());
						fprintf(fisicv,"%f ",stdev / mean); // ISI coefficient of variation for that neuron
					}
					
					// per pop
					fprintf(fspikes_pop[pop_belonging],"%i ",spikes);
					if(neuron_isi[ii].size() > 1) {
						double total = std::accumulate(neuron_isi[ii].begin(),neuron_isi[ii].end(),0.0);
						double mean = total / neuron_isi[ii].size();
						std::vector<double> diff(neuron_isi[ii].size());
						std::transform(neuron_isi[ii].begin(), neuron_isi[ii].end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
						double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
						double stdev = sqrt(sq_sum / neuron_isi[ii].size());
						fprintf(fisicv_pop[pop_belonging],"%f ",stdev / mean); // ISI coefficient of variation for that neuron
					}
					
				}
				
				for(int ii=0; ii < pops.size(); ii++) {
					printf("Population %s average firing rate: %f\n",pops[ii]->name.c_str(),((float)spikes_pop[ii] * 1000.0f / t_end)/(float)pops[ii]->num_neurons);
					fclose(fintervals_pop[ii]);
					fclose(fspikes_pop[ii]);
					fclose(fisicv_pop[ii]);
				}
				free(fisicv_pop);
				free(fintervals_pop);
				free(fspikes_pop);
				fclose(fspikes);
				fclose(fintervals);
				fclose(fisicv);
				fclose(fneuron_activity);
		#endif
			}
			
			
			
			if(globalNeuronRefs != NULL) free(globalNeuronRefs);
			
			delete communicator;
			delete partition;
			
			
			if(custom_partitioning != NULL) free(custom_partitioning);

			free(computation_times);
			free(propagate_times);
			free(sync_time);
			free(idle_time);
		#ifdef RECORD_PROCESS_TRACE
			free(trace_timings);
		#endif
			free(comp_times);
			free(comm_sent);
			free(p_times);
			free(sim_times);
			free(idle_times);
			free(s_times);
			free(comp_variances);
			free(idle_variances);
			free(sync_variances);
			free(runtime_neighbours);
			free(non_empty_comm_steps);

		}

	}

	// clean up
	for(int ii = 0; ii < pops.size(); ii++) {
		delete pops[ii];
	}
	
	for(int ii=0; ii < injections.size(); ii++) {
		delete injections[ii];
	}
	for(int ii=0; ii < cell_params.size(); ii++) {
		free(cell_params[ii]);
	}
	
	MPI_Finalize();
	return 0;
}
