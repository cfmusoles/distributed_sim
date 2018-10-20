#ifndef STATIC_CONNECTIVITY__H
#define STATIC_CONNECTIVITY__H

#include "Connectivity.h"
#include "Utils.h"
#include <algorithm>

class StaticConnectivity : public Connectivity {
public:
	
	StaticConnectivity(Population* p) : Connectivity(p,p) {
	}
	virtual ~StaticConnectivity() {}
	
	virtual void generate_connectivity(int pop_size) {
		// connectivity is statically defined elsewhere (by loading file)
		// TODO centralise process in this class
		return;
		
	}
	
	//virtual bool has_connection(int from, int to) {
	//	return std::find(connections[from].begin(),connections[from].end(),to) != connections[from].end();
	//}
	
};

#endif
