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
	Graph *graph = new Graph(argv[1]);

	bool print = false;
	if(strcmp(argv[argc-1], "output") == 0) print = true;

	if(argc > 2&&strcmp(argv[2], "list-heap") == 0) {
		graph->core_decomposition(1, print);
	}
	else if(argc > 2&&strcmp(argv[2], "array-heap") == 0) {
		graph->core_decomposition(2, print);
	}
	else if(argc > 2&&strcmp(argv[2], "h-index") == 0) {
		graph->core_decomposition_hindex(print);
	}
	else if(argc > 2&&strcmp(argv[2], "hierarchy") == 0) {
		graph->core_hierarchy(print);
	}
	else graph->core_decomposition(0, print);

	printf("**********************************\n");

	return 0;
}

void print_usage() {
	printf("Usage: [1]exe [2]graph-dir [3]algorithm [4 optional] output\n");
	printf("    Algorithms: list-heap, array-heap, h-index, null (i.e., nothing), hierarchy\n");
}
