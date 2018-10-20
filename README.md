# Distributed SNN simulation

This repository contains code associated with the paper submitted to Frontiers in Neuroinformatics titled ***"Communication sparsity in distributed Spiking Neural Network Simulations to improve scalability"*** (approval pending). It consists of a simple SNN timestepped simulator and the implementation of two SNN models: Cortical Microcircuit (Potjans and Diesmann, 2014), and  Macaque Visual Cortex (Schmidt, 2018). 

Both models are implemented in C++, where different neuron allocation and communication strategies are evaluated. The main features of the implementation are:

* Time-step driven simulation.
* Single neuron (leak integrate and fire) and synaptic models (exponential). 
* Distributed across computing nodes (processes). Within a process, computation is performed sequentially (computation performance is out of the scope of this work).
* MPI library used for interprocess communication.

These features are built to be representative of the wider class of SNN simulators, with functionally equivalent output (network activity) for the same model execution. Thus, the findings on this work could inform development in other time-step, distributed SNN simulators (such as NEST, NEURON, Arbor).

Profiling timing is done internally wrapping specific functions with MPI_Wtime calls. We are mainly interested in measuring time spent in computation, in implicit synchronisation and in data exchange, as well as how these impact overall performance.

* Computation: includes the sequential computation of neuron and synapse state updates, as well as effecting the received spiking information from other processes to local neurons.
* Implicit synchronisation: to quantify load balance, a global artificial barrier is introduced at the end of the computation cycle to measure waiting time on each process\footnote{This is implemented with MPI\_Barrier calls. Although there are limitations to this approach (a process is only held until it is notify that all other processes have entered the barrier; but message propagation delays may result in additional imbalance) it is a good approximation of time to synchronise.}. 
* Data exchange: measurement of the time it takes for the selected communication pattern to initiate and complete spike propagation.
* Simulation time: overall global time to complete the simulation, including computation, implicit synchronisation and data exchange.


## Results in the paper

Results files from experiments in paper stored in results/archer and results/azure
Python scripts to generate the figures in paper in plots/
Implementation of communication algorithms in paper
* PEX --> include/PEXCommunicator.h
* NBX --> include/NBXCommunicator.h
Implementation of workload distribution algorithms
* Random Balanced --> include/RandomBalancedPartitioning.h
* Round Robin --> include/RoundRobinPartitioning.h
* Hypergraph partitioning --> include/HypergraphPartitioning.h
- Description of main files
* src/distSim.cpp --> simulator code, model builder
Command used to get results in the paper

For the scaling experiment of the CM model
> mpirun -n $NUM_PROC -ppn $PROCESS_PER_NODE distSim -n cm_test -c $COMM_PATTERN -p $PARTITIONING -s $SEED -k 500 -f 1000 -t 750 -m "cm"

For the scaling experiment of the MVC model
> mpirun -n $NUM_PROC -ppn $PROCESS_PER_NODE distSim -n cm_test -c $COMM_PATTERN -p $PARTITIONING -s $SEED -k 1000 -f 200 -t 350 -m "mvc"

Where (explain the variable names)


## Instructions on how to run

### Dependencies

The C++ project dependencies are: 
* [metis 5.1](http://glaros.dtc.umn.edu/gkhome/metis/metis/download)
* [parmetis 4.0.3]( http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)
* [zoltan 3.83](http://www.cs.sandia.gov/Zoltan/Zoltan_download.html)
* An instalation of MPI (OpenMPI 3.0.0 or MPI Intel 5.X used in the paper)
* C++11 minimum

### Compilation 

The project can be compiled using the `Makefile` included. First, modify the `Makefile` path variables to dependency libraries to point to the appropriate folder in your instalation:
```
METIS_PATH=$(HOME)/metis_build
PARMETIS_PATH=$(HOME)/parmetis_ompi_build
ZOLTAN_PATH=$(HOME)/Zoltan_v3.83/build
```
Then you can compile the project with
> make

	
### Running. Parameters
```
-n: test name
-c: communication pattern ('pex' or 'nbx')
-p: workload distribution (partitioning of neurons) ('randomBalanced','roundrobin' or 'hypergraphPartitioning')
-s: random seed
-w: neuronal activity file 
-k: 0 to 1000 (1000 = 100%) fraction of synapses (scaling model)
-f: 0 to 1000 (1000 = 100%) fraction of neurons (scaling model)
-t: simulated time in ms
-m: model selection ('cm' or 'mcv')
```

## References
- Schmidt, M., Bakker, R., Hilgetag, C. C., Diesmann, M., & van Albada, S. J. (2018). Multi-scale account of the network structure of macaque visual cortex. Brain Structure and Function, 223(3), 1409–1435. https://doi.org/10.1007/s00429-017-1554-4
- Potjans, T. C., & Diesmann, M. (2014). The cell-type specific cortical microcircuit: Relating structure and activity in a full-scale spiking network model. Cerebral Cortex, 24(3), 785–806. https://doi.org/10.1093/cercor/bhs358

