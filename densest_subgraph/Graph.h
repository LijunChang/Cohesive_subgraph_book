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

using lint = long long;
#define min(a,b) ((a)<(b)?(a):(b))

const lint LINT_INF = (lint)2000000000*(lint)2000000000;

class Graph {
private:
	std::string dir; //input graph directory
	ui n; //#nodes of the graph
	ui m; //#edges of the graph

	ui original_n; //#nodes of the input graph
	ui *component_id; //component ids of SCCs of nodes

	ui *pstart; //start positions of neighbors of vertices in the array "edges"
	ui *pend; //end positions of neighbors of vertices in the array "edges"
	ui *edges; //concatenation of neighbors of all vertices

	ui *revers; // linke two directions of an undirected edge together
	lint *capacity_pre; // initial capacities of edges
	lint *capacity; // capacities of edges
	ui *level; // used in dinic
	ui *now; // used in dinic
	ui *que; // used in dinic

	lint total_flow_pre, lambda_pre;

	std::vector<ui> res;

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph_binary() ;
	void read_graph(const char *input_file) ;

	// compute the maximal densest subgraph
	void densest_subgraph() ;

	void print_densest_subgraph() ;
	void print_densest_subgraph(const char *output_file) ;

private:
	// core decomposition without invoking data structures
	ui core_decomposition(ui *peel_sequence, ui *core, const ui n, const ui *pstart, const ui *pend, const ui *edges) ;

	void compute_SCCs(ui &scc_n, ui *scc_starts, ui *scc_ids, const ui n, const ui *pstart, const ui *pend, const ui *edges) ;

	void construct_lambda_graph() ;

	void greedy_and_reduce(ui &best_n, ui &best_m, ui *peel_sequence, ui *core, std::vector<ui> &ids, ui *degree, ui *rid) ;

	bool bfs() ;
	lint dinic(ui u, lint flow) ;
	lint min_cut(char initialize, lint lambda) ;
	void save_S(std::vector<ui> &res, char *vis) ;
};

#endif /* GRAPH_H_ */
