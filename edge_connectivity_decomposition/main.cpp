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
	if(argc < 3) {
		print_usage();
		return 0;
	}
	Graph *graph = new Graph(argv[1]);

	bool print = false;
	if(strcmp(argv[argc-1], "output") == 0) print = true;

	if(strcmp(argv[2], "kecc") == 0) {
		if(argc < 4||(argc == 4&&print)) print_usage();
		else graph->k_edge_connected_component((ui)atoi(argv[3]), print);
	}
	else if(strcmp(argv[2], "kecc-space") == 0) {
		if(argc < 4||(argc == 4&&print)) print_usage();
		else graph->k_edge_connected_component_space((ui)atoi(argv[3]), print);
	}
	else print_usage();

	printf("**********************************\n");

	return 0;
}

void print_usage() {
	printf("Usage: [1]exe [2]graph-dir [3]\"kecc\" or \"kecc-space\" [4]k [5 optional]\"output\"\n");
}
