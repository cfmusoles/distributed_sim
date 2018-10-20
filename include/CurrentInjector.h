#ifndef CURRENTINJECTOR__H
#define CURRENTINJECTOR__H

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "Injector.h"


class CurrentInjector : public Injector {
public:
	float fixed_current; // fixed injection current per time step									
	
	CurrentInjector(int sz, float fc, Population* pop, float s, float e, float p, float dt, float k_scale ) {
		fixed_current = fc;
		init(pop, p, sz, dt, k_scale, e, s);
	}
	virtual ~CurrentInjector() {}
	
protected:
	virtual float calculateCurrent(float current_time, float last_spike) {
		return fixed_current;
	}
	
	virtual float calculateNextSpike(float current_time) {
		return current_time; // force injection on every time step
	}
};

#endif
