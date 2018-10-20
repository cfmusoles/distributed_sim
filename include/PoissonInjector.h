#ifndef POISSON_INJECTOR__H
#define POISSON_INJECTOR__H

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "Injector.h"


class PoissonInjector : public Injector {
public:
	double frequency;									// spike frequency (1 / times per millisecond) 8 Hz = 1 / 0.008
	
	PoissonInjector(int sz, int freq, Population* pop, float s, float e, float p, float dt, float k_scale ) {
		frequency = 1.0f / (1000.0f / freq);
		init(pop, p, sz, dt, k_scale, e, s);
	}
	virtual ~PoissonInjector() {}
	
protected:
	virtual float calculateCurrent(float current_time, float last_spike) {
		return  w * expf(-(current_time-last_spike) / TAU_S );
	}
	
	virtual float calculateNextSpike(float current_time) {
		// Poisson cummulative distribution http://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/
		return current_time + -logf(1.0f - (double) random() / (double)RAND_MAX) / frequency;
	}
};



#endif
