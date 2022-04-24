/*
 * main.cpp
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#include "Graph.h"

void print_usage() ;

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("!!! You may want to define NDEBUG in utilities/Defines.h to get better performance!\n");
#endif
	if(argc < 2) {
		print_usage();
		return 0;
	}
#ifndef NDEBUG
	printf("**** DensestSubgraph (Debug) build at %s %s ***\n", __TIME__, __DATE__);
#else
	printf("**** DensestSubgraph (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

	Graph *graph = new Graph(argv[1]);

	graph->read_graph_binary();
	//graph->read_graph(argv[1]);

	graph->densest_subgraph();

	//graph->print_densest_subgraph(argv[2]);
	if(strcmp(argv[argc-1], "output") == 0) graph->print_densest_subgraph();

	printf("**********************************\n");

	return 0;
}

void print_usage() {
	printf("Usage: [1]exe [2]graph-dir [3 optional] \"output\"\n");
}
