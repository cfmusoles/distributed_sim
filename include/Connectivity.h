#ifndef CONNECTIVITY__H
#define CONNECTIVITY__H

#include <vector>

class Connectivity {
public:
	Population* from;
	Population* to;
	std::vector<std::vector<int> > connections;
	
	Connectivity(Population* f, Population* t) {
		from = f;
		to = t;
		connections.resize(from->num_neurons);
	}
	virtual ~Connectivity() {}
	
	virtual void generate_connectivity(int pop_size) { }
	
	//virtual bool has_connection(int from, int to) {
	//	return false;
	//}
	
	virtual void clean_data_structs() {
		connections.clear();
		connections.swap(connections);
	}
};

#endif
