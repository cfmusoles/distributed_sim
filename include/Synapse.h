#ifndef SYNAPSE__H
#define SYNAPSE__H

#include <stdint.h>
#include <math.h>
#include "Neuron.h"
#include "Utils.h"
#include "NeuronParams.h"


class Synapse {
public:
	//int pre_neuron;
	double t_f;					// last time fired
	double w;					// synaptic strength or weight
	Neuron* post_syn; 			// postsynaptic neuron pointer
	uint32_t spike_buffer;		// spike buffer (to allow tracking of multiple spikes and delay)
	short delay_ticks;			// number of sim ticks that represent the synaptic delay
	double tau_s;
	
	Synapse() {}
	Synapse(int pren, Neuron* postsynaptic, bool isExcitatory, float dt, NeuronParams* c, float k_scale, float weight = 0.0f, float delay = 0.0f) {
		init(pren,postsynaptic, isExcitatory, dt, c, k_scale,weight,delay);
	} 
	virtual ~Synapse() {}
	
	void init(int pren, Neuron* postsynaptic, bool isExcitatory, float dt, NeuronParams* c, float k_scale, float weight = 0.0f, float delay = 0.0f) {
		//pre_neuron = pren;
		post_syn = postsynaptic;
		init_synapse(isExcitatory,dt,k_scale,c,weight,delay);
	}
	
	void spike(int time_steps_ago) {
		spike_buffer = spike_buffer | (1 << (delay_ticks - time_steps_ago));
	}
	
	void update_current(float current_t) {
		if(spike_buffer & 0x00000001) {
			// deliver spike
			t_f = current_t;
		}
		post_syn->i_e += w * expf(-(current_t-t_f) / tau_s);
		spike_buffer = spike_buffer >> 1;
	}
	
	bool operator < (const Synapse& syn2) const
    {
        return (post_syn->id < syn2.post_syn->id);
    }

private:
	void init_synapse(bool isExcitatory, float dt, float k_scale, NeuronParams* cell_params, float weight, float delay) {
		/*w = randgauss(cell_params->w_exc_mean, cell_params->w_exc_dev);
		w = w < 0 ? 0 : w;
		float d = dt;
		if(isExcitatory) {
			d = randgauss(cell_params->exc_delay_mean, cell_params->exc_delay_dev);
			tau_s = cell_params->exc_tau_s;
		} else {
			w *= cell_params->w_inh_g;
			d = randgauss(cell_params->inh_delay_mean, cell_params->inh_delay_dev);
			tau_s = cell_params->inh_tau_s;
		}
		d = d < 0 ? 0 : d;
		d -= remainderf(d,dt); // round off delay to multiple of dt
		delay_ticks = d/dt;
		w /= sqrt(k_scale);
		t_f = -10;
		spike_buffer = 0x0L;*/
		
		w = weight;
		if(isExcitatory) {
			tau_s = cell_params->exc_tau_s;
		} else {
			tau_s = cell_params->inh_tau_s;
		}
		double d = delay < 0 ? 0 : delay;
		d -= remainderf(d,dt); // round off delay to multiple of dt
							
		delay_ticks = d/dt;
		w /= sqrt(k_scale);
		t_f = -10;
		spike_buffer = 0x0L;
		
		
	}
};
#endif
