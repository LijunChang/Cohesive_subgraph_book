/*
 * main.cpp: transform a graph from edge list to our binary represenation
 *
 *  Created on: 21Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#include "../utilities/Defines.h"
#include "../utilities/Utility.h"

using namespace std;

int main(int argc, char *argv[]) {
	if(argc < 3) {
		printf("Usage: [1]exe [2]graph-dir [3]edge-list-filename\n");
		printf("\tIn the edgelist file, lines starting with non-number characters are skipped!\n");
		return 0;
	}

	vector<ui> nodes;
	vector<pair<ui,ui> > edges;

	string dir = string(argv[1]);
	string file_name = string(argv[2]);
	FILE *f = Utility::open_file((dir + string("/") + file_name).c_str(), "r");

	char buf[1024];
	ui a, b;
	while(fgets(buf, 1024, f)) {
		char comment = 1;
		for(ui j = 0;buf[j] != '\0';j ++) if(buf[j] != ' '&&buf[j] != '\t') {
			if(buf[j] >= '0'&&buf[j] <= '9') comment = 0;
			break;
		}
		if(comment) continue;

		for(ui j = 0;buf[j] != '\0';j ++) if(buf[j] < '0'||buf[j] > '9') buf[j] = ' ';
		sscanf(buf, "%u%u", &a, &b);
		if(a == b) continue;
		nodes.push_back(a);
		nodes.push_back(b);
		edges.push_back(make_pair(a,b));
		edges.push_back(make_pair(b,a));
	}

	fclose(f);

	sort(nodes.begin(), nodes.end());
	nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

	printf("min id = %u, max id = %u, n = %lu\n", nodes.front(), nodes.back(), nodes.size());

	sort(edges.begin(), edges.end());
	edges.erase(unique(edges.begin(), edges.end()), edges.end());

	map<ui,ui> M;
	for(ui i = 0;i < nodes.size();i ++) M[nodes[i]] = i;

	char preserved = 1;
	for(ui i = 0;i < nodes.size();i ++) if(nodes[i] != i) preserved = 0;
	if(!preserved) printf("Node ids are not preserved!\n");

	ui n = nodes.size();
	ui m = edges.size();
	printf("n = %s, (undirected) m = %s\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *pstart = new ui[n];
	ui *edge = new ui[m];

	ui j = 0;
	for(ui i = 0;i < n;i ++) {
		pstart[i] = j;
		while(j < m&&edges[j].first == nodes[i]) {
			edge[j] = M[edges[j].second];
			++ j;
		}
	}
	pstart[n] = j;

	f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "wb");

	ui tt = sizeof(ui);
	fwrite(&tt, sizeof(ui), 1, f);
	fwrite(&n, sizeof(ui), 1, f);
	fwrite(&m, sizeof(ui), 1, f);

	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1]-pstart[i];
	fwrite(degree, sizeof(ui), n, f);
	fclose(f);

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "wb");
	fwrite(edge, sizeof(ui), m, f);
	fclose(f);

	delete[] pstart;
	delete[] edge;
	delete[] degree;
}
