#ifndef UTILS__H
#define UTILS__H

#include <stdlib.h>
#include <math.h>

#include <cstdio>

#ifdef VERBOSE
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...) 
#endif

#define PI 3.141592654

template<typename T>
int convert_to_bytes(T value, char* buffer) {
	for (int ii = 0; ii < sizeof(T); ++ii)
	{
		buffer[ii] = (value >> 8 * ii) & 0xFF;
		//printf("%02X\n",(unsigned char)buffer[ii]);
	}
	return sizeof(T);
}

template<typename T>
int convert_from_bytes(char* buffer, T* value) {
	*value = 0;
	for (int ii = 0; ii < sizeof(T); ++ii)
	{
		unsigned long partial_value = (unsigned long)((unsigned char)buffer[ii]) << (8 * ii);
		//printf("Partial: %u\n",partial_value);
		*value |= partial_value;
	}
	return sizeof(T);
	//printf("Reconstructed: %lu\n",*value);
}

double randgauss(float mean, float std, bool reset = false)
{
	static double U, V;
	static int phase = 0;
	double Z;

	if(reset) {
		phase = 0;
		return -1;
	} 
	
	if(std > 0) {
		
		if(phase == 0) {
			U = (rand() + 1.) / (RAND_MAX + 2.);
			V = rand() / (RAND_MAX + 1.);
			Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
		} else
			Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

		phase = 1 - phase;

		return Z*std + mean;
	} else {
		return mean;
	}
}

double rand01() {
	return (double)rand() / (double)RAND_MAX;
}

double clamp(double value, double min, double max) {
	return value < max ? ( value > min ? value : min) : max;
}

double max(double value1, double value2) {
	return value1 > value2 ? value1 : value2;
}

double min(double value1, double value2) {
	return value1 < value2 ? value1 : value2;
}

unsigned int encodeSpike(int spike_timing, int neuronId) {
	unsigned int message = neuronId | (spike_timing << 27);
	return message;
}

void decodeSpike(unsigned int message, int *spike_timing, int* neuronId) {
	*neuronId = 0x07ffffff & message;
	*spike_timing = message >> (27);
}

#endif
