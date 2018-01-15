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
#include "../data_structures/UnionFind.h"

struct Edge {
	Edge *pre, *next;
	Edge *reverse;
	ui vertex;
};

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

	// compute k-edge connected components
	void k_edge_connected_component(ui K, bool print) ;
	// compute k-edge connected components in a space-effective manner
	void k_edge_connected_component_space(ui K, bool print) ;

private:
	void read_graph_binary() ;

	// k-core based pruning
	void k_core_prune(ui K, ui *Q, ui Q_n, char *computed, ui *degree, ui *pend) ;
	void print_kecc(ui K, ui c_n, ui *cstart, ui *ids) ;

	// initialize the partition graph
	void initialize_pgraph(ui s,ui *Q, char *vis, char *computed, ui *pend, ui *pend_local, UnionFind *UF, ui *representative, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last) ;
	ui decomposition(const ui s, const ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, ui *pend_local, ListLinearHeap *heap, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last, UnionFind *UF, ui *representative) ;
	void remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *pend, ui *Q, char *computed, ui *degree, UnionFind *UF) ;

	// construct the partition graph
	void construct_pgraph(ui s, ui *Q, char *vis, char *computed, ui *pend, Edge **graph_head, Edge *graph_edges, ui *sv_next, ui *sv_last) ;
	ui decomposition(ui s, ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, Edge **graph_head, Edge *graph_edges, ListLinearHeap *heap, ui *sv_next, ui *sv_last) ;
	void merge(Edge **graph_head, Edge *edges, ui u, ui v, ui *keys, ui *sv_next, ui *sv_last) ;
	void add_edge(ui u, Edge *e, Edge **graph_head) ;
	void delete_edge(ui u, Edge *e, Edge **graph_head) ;
	void remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *component_id, ui *pend, ui *Q, char *computed, ui *degree) ;
};

#endif /* GRAPH_H_ */
