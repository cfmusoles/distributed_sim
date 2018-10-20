#ifndef NEURON__H
#define NEURON__H

#include <math.h>
#include <Utils.h>
#include "NeuronParams.h"
#include <vector>

// 1 F = 10^12 pF
// 1 MOhm = 1000000 Ohm
// TAU_M = R * Cm -> SECONDS = ohms * farads -> 0.01 s = Rm(ohms) * 250 * 10^-12 F -> Rm = 4*10^7 ohms = 40

class Neuron {
public:
	int id;					// global ID
	double v;				// membrane voltage or potential
	double i_e;				// input current injected (summation) 
	NeuronParams* cell_params;
	
	Neuron() {}
	
	Neuron(int gid, NeuronParams* c) {
		init(gid,c,0);
	}
	virtual ~Neuron() {}
	
	void init(int gid, NeuronParams* c, float init_v) {
		id = gid;
		cell_params = c;
		v = init_v;
		last_fired = -cell_params->tau_ref;
		i_e = 0;
	}
	
	bool update_potential(float dt, float current_t) {  
		bool fired = false;
		if(current_t >= cell_params->tau_ref + last_fired) { // refractory period after firing
			// update membrane voltage (LIF)
			i_e += cell_params->i_offset;
			// from matlabhttp://neuroscience.ucdavis.edu/goldman/Tutorials_files/Integrate%26Fire.pdf
			double v_inf = cell_params->v_rest + i_e * cell_params->r_m;
			v = v_inf + (v - v_inf) * expf(-dt/cell_params->tau_m);
			
			// check for spikes
			if(v > cell_params->v_thres) {
				v = cell_params->v_reset;
				last_fired = current_t;
				fired = true;
			}
		}
		// clear synaptic current
		i_e = 0;
		return fired;
	}
	
private:
	double last_fired;
	
};

#endif
