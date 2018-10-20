#ifndef INJECTOR__H
#define INJECTOR__H

#include <algorithm>
#include "Simulation.h"

#define DELAY 1.5f		// synaptic delay for current injection spikes
#define TAU_S 0.5f
#define W_EXC 0.0878f		// in nA


class Injector {
public:
	Population* population;
	float p_conn;
	
	int size;
	float start;
	float end;
	short delay_ticks;									// number of sim ticks that represent the synaptic delay
	float w;											// synaptic strength or weight
	std::vector<double> last_spike;
	std::vector<double> next_spike;
	std::vector<uint32_t> spike_buffer;					// spike buffer (to allow tracking of multiple spikes and delay)
	std::vector<float> i_e;								// current to be injected at current time
	std::vector<std::vector<int> > post_syn_neuron_idx;	// list of postsyn neurons targeted by each current injector "neuron" id 
	
	Injector() { }
	virtual ~Injector() {}
	
	void create_connections(int process_id) {
		// create injection connections for the entire network
		for(int jj=0; jj < size; jj++) {
			for(int kk=0; kk < population->neuron_ids.size(); kk++) {
				if(rand01() <= p_conn) {// probabilistic connectivity (not all to all)
					post_syn_neuron_idx[jj].push_back(population->neuron_ids[kk]);
				}
			}
			/*for(int kk=population->first_idx; kk < population->first_idx + population->num_neurons; kk++) {
				if(rand01() <= p_conn) {// probabilistic connectivity (not all to all)
					post_syn_neuron_idx[jj].push_back(kk);
				}
			}*/
		}
	}
	
	void cull_to_local_post_syn_neurons(idx_t* partitioning, int process_id) {
		// remove any connections to non local neurons
		for(int jj=0; jj < size; jj++) {
			std::vector<int> vec;
			for(int ii=0; ii < post_syn_neuron_idx[jj].size(); ii++) {
				if(partitioning[post_syn_neuron_idx[jj][ii]] == process_id) {
					vec.push_back(post_syn_neuron_idx[jj][ii]);
				}
			}
			std::swap(vec,post_syn_neuron_idx[jj]);
			
		}
	}
	
	void update(double current_time) {
		// loop through population reservoire
		for(int ii=0; ii < size; ii++) {
			if(current_time >= start && current_time <= end) {
				// check spike buffer
				if(spike_buffer[ii] & 0x01) {
					// record spike
					last_spike[ii] = current_time;
				}
				// update spike buffer
				if(next_spike[ii] <= current_time) {
					if(next_spike[ii] > 0) {
						spike_buffer[ii] = spike_buffer[ii] | (1 << delay_ticks);
					}
					next_spike[ii] = calculateNextSpike(current_time);
				}
			}
			i_e[ii] = calculateCurrent(current_time,last_spike[ii]);
			spike_buffer[ii] = spike_buffer[ii] >> 1;
			//if(i_e[ii] > 0.00001) printf("Next spike: %f, last spike: %f, spike_buffer: %i; current: %f\n",next_spike[ii],last_spike[ii],spike_buffer[ii], i_e[ii]);
		}
	}
protected:
	
	void init(Population* pop, float p, int sz, float dt, float k_scale, float e, float s) {
		delay_ticks = DELAY/dt;
		w = W_EXC / sqrt(k_scale);
		start = s;
		end = e;
		population = pop;
		p_conn = p;
		
		size = sz;
		next_spike.resize(size);//-1.0f;
		last_spike.resize(size);// = -10;
		spike_buffer.resize(size);// = 0x0L;
		i_e.resize(size);
		post_syn_neuron_idx.resize(size);
		for(int ii=0; ii < size; ii++) {
			next_spike[ii] = -1.0f;
			last_spike[ii] = -10;
			spike_buffer[ii] = 0x0L;
			i_e[ii] = 0.0f;
		}
	}
	
	virtual float calculateCurrent(float current_time, float last_spike) { return 0.0f;}
	virtual float calculateNextSpike(float current_time) { return 0.0f; }
};

#endif
