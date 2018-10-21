# Description of result files

The files contained in this folder correspond to results generated in the [ARCHER supercomputer](http://archer.ac.uk/about-archer/) for the Cortical Microcircuit scaling experiments.

Each experiment contains two types of result files:
* Result timings: contain global simulation metrics, such as simulation time, synchronisation time, number of spikes, average number of runtime neighbours, etc. 
* Tracing files: contain timings per process, such as computation, synchronisation, idle time, etc.

## Result timing files

The result timings files are named as follows:
> sparse_comm_st_$DISTRIBUTION_$COMM_PATTERN__$TOTAL_NUM_PROC

The result timing files contain a header that describes the values represented on each column. The more representative colums are:
* Build time: time it took to build the model and prepare the simulation (workload distribution)
* Sim time: time it took to complete the simulation
* Propagation time: time the simulation spent in communication of spikes (including data exchange and implicit synchronisation)
* Sync time: time the simulation spent in data exchange
* Idle time: time the simulation spent in implicit synchronisation
* Bytes sent: minimum bytes exchanged between processes (does not include overhead of MPI communication)
* Edgecut: percentage of synapses that connect neurons hosted in two different processes
* Spikes sent: number of spikes copies sent across processes (one spiking neuron needs to send as many copies as target processes it connects to). 
* Remote spikes: number of spikes that had to be copied during simulation
* Average runtime neighbours: on average, number of processes each process communicates with during each communication phase.

## Tracing files

For each scaling experiment, there are as many tracing files as number of processes used. The files are named as follows:
> sparse_comm_st_$DISTRIBUTION_$COMM_PATTERN_$TOTAL_NUM_PROC_$PROCESS_ID

The tracing files contain a header that describes the values represented on each column. The more representative colums are:
* Comp time: time the process spent in computation phase
* Idle time: time the process spent in implicit synchronisation phase
* Sync time: time the process spent in data exchange phase
* Runtime neighbours: average number of processes the process communicates with during each communication phase.
