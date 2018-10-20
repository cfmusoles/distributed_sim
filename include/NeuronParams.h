#ifndef NEURON_PARAMS_H
#define NEURON_PARAMS_H

struct NeuronParams {
	// neuron params
	float v_start_mean;		// mean start voltage (mV)
	float v_start_dev;		// std for start voltage (mV)
	float tau_m;			// membrane time constant (ms) 
	float r_m;				// membrane resistance (MOhm) 
	float tau_ref;			// refractory period after spike (ms) 
	float v_rest;			// resting potential (mV)
	float v_thres;			// spike threshold (mV) 
	float v_reset;			// voltage after spike (mV) 
	float v_spike;			// nominal spiking potential (mV) (only for drawing purposes) 
	float i_offset;			// constant current input (nA)
	float c_m;				// cm (nF)
	// synapse params
	float w_exc_mean;		// nA (PyNN default is nA)
	float w_exc_dev;			// nA (PyNN default is nA)
	float w_inh_g;			// g is a multiplier factor with respect to W_EXC
	float exc_tau_s;			// time constant for exc synapses (ms) 
	float inh_tau_s;		// time constant for inh synapses (ms) 
	float exc_delay_mean;	// delay on exc synaptic propagation (ms)
	float exc_delay_dev;		// delay on exc synaptic propagation (ms)
	float inh_delay_mean;	// delay on inh synaptic propagation (ms)
	float inh_delay_dev;		// delay on inh synaptic propagation (ms)
};

#endif
