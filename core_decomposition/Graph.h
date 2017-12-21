/*
 * Graph.h
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include "../utilities/Defines.h"
#include "../utilities/Utility.h"
#include "../utilities/Timer.h"
#include "../data_structures/ListLinearHeap.h"
#include "../data_structures/ArrayLinearHeap.h"
#include "../data_structures/UnionFind.h"

class Graph {
private:
	std::string dir; //input graph directory
	ui n; //#nodes of the graph
	ui m; //#edges of the graph

	ui *pstart; //start positions of neighbors of vertices in the array "edges"
	ui *edges; //concatenation of neighbors of all vertices

public:
	Graph(const char *_dir) ;
	~Graph() ;

	// linear heap based core decomposition
	void core_decomposition(ui alg, bool print) ;
	// hindex based core decomposition
	void core_decomposition_hindex(bool print) ;

	void core_hierarchy(bool print) ;

private:
	void read_graph_binary() ;

	// for each distinct key k, compute the size of the gcc for the graph induced by vertices whose keys are no less than k
	void compute_core_gccs(const ui *seq, const ui *key) ;

	// core decomposition without invoking data structures
	ui core_decomposition(ui *peel_sequence, ui *core) ;
	// core decomposition with listlinearheap or arraylinearheap
	ui core_decomposition_linear_heap(ui alg, ui *peel_sequence, ui *core) ;
	// core decomposition based on hindex
	ui core_decomposition_hindex(ui *core) ;
	void print_cores(const ui *core) ;
};

#endif /* GRAPH_H_ */
