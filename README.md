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

The application must be run with a series of arguments that configure the simulation:
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
For example, the following command would be used to run a 100ms simulation of the full CM model, with NBX communication algorithm and round robin distribution, in 12 distributed processes, with 123456 as random seed:
> mpirun -n 12 distSim -n test -c nbx -p roundrobin -t 100 -m cm -s 123456

## Results in the paper

The result files presented in the paper for the CM simulations are included in the folder `results/archer`, and the MVC simulations are included in `results/azure`.

All the figures were generated with Python scripts, which are included in the folder `plots`.

The code implementation of the different communication algorithms presented in the paper can be found in:
* PEX --> `include/PEXCommunicator.h`
* NBX --> `include/NBXCommunicator.h`

The code implementation of the different workload distribution algorithms presented in the paper can be found in:
* Random Balanced --> `include/RandomBalancedPartitioning.h`
* Round Robin --> `include/RoundRobinPartitioning.h`
* Hypergraph partitioning --> `include/HypergraphPartitioning.h`

The code for the simulator and the model builder can be found in `src/distSim.cpp`

### Command used to get results in the paper

The following command runs the scaling experiment of the CM model:
> mpirun -n $NUM_PROC -ppn $PROCESS_PER_NODE distSim -n cm_test -c $COMM_PATTERN -p $PARTITIONING -s $SEED -k 500 -f 1000 -t 750 -m "cm"

The following command runs the scaling experiment of the MVC model:
> mpirun -n $NUM_PROC -ppn $PROCESS_PER_NODE distSim -n cm_test -c $COMM_PATTERN -p $PARTITIONING -s $SEED -k 1000 -f 200 -t 350 -m "mvc"

In both cases, the variable names (in capitals and preceded by the symbol $) are:
* `$NUM_PROC`: number of distributed (MPI) processes used in the simulation (96,192,384,768,1536,3072)
* `$PROCESS_PER_NODE`: number of processes in a computing node (24 for ARCHER, 16 for Azure H16mr nodes)
* `$COMM_PATTERN`: the communication pattern used in the simulation (`'randomBalanced'`, `'roundrobin'` or `'hypergraphPartitioning'`
* `$SEED`: random seed of the experiment (random numbers used for repeatability)

The values used to produce the results in the paper are indicated in brackets.

## Instructions on how to run

### Dependencies

The C++ project dependencies are: 
* [metis 5.1](http://glaros.dtc.umn.edu/gkhome/metis/metis/download)
* [parmetis 4.0.3]( http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)
* [zoltan 3.83](http://www.cs.sandia.gov/Zoltan/Zoltan_download.html)
* An instalation of MPI (OpenMPI 3.0.0 or MPI Intel 5.X used in the paper)
* C++11 minimum


## References
- Schmidt, M., Bakker, R., Hilgetag, C. C., Diesmann, M., & van Albada, S. J. (2018). Multi-scale account of the network structure of macaque visual cortex. Brain Structure and Function, 223(3), 1409–1435. https://doi.org/10.1007/s00429-017-1554-4
- Potjans, T. C., & Diesmann, M. (2014). The cell-type specific cortical microcircuit: Relating structure and activity in a full-scale spiking network model. Cerebral Cortex, 24(3), 785–806. https://doi.org/10.1093/cercor/bhs358

