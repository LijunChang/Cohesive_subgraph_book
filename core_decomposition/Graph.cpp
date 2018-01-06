/*
 * Graph.cpp
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#include "Graph.h"

using namespace std;

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = 0;
	m = 0;

	pstart = nullptr;
	edges = nullptr;
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
}

void Graph::core_decomposition(ui alg, bool print) {
	if(alg > 2) {
		printf("Parameter [alg] in core_decomposition is wrong!\n");
		return ;
	}

	char name[3][16] = {"", "list-heap", "array-heap"};
#ifdef NDEBUG
	printf("*** core_decomposition %s (Release): %s ***\n", name[alg], dir.c_str());
#else
	printf("*** core_decomposition %s (Debug): %s ***\n", name[alg], dir.c_str());
#endif

	read_graph_binary();

	Timer timer;
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n]; //core[i] is the core value of vertex i

	ui max_core = 0;
	if(alg == 0) max_core = core_decomposition(peel_sequence, core);
	else max_core = core_decomposition_linear_heap(alg, peel_sequence, core);

	//compute_core_gccs(peel_sequence, core);

	ui two_core = n;
	for(ui i = 0;i < n&&core[peel_sequence[i]] < 2;i ++) -- two_core;
	printf("\tMax core value: %u, two_core size: %s\n", max_core, Utility::integer_to_string(two_core).c_str());

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) if(i == 0||core[peel_sequence[i]] != core[peel_sequence[i-1]]) {
		printf("Core: %u, #vertices: %s\n", core[peel_sequence[i]], Utility::integer_to_string(n - i).c_str());
	}
#endif

	if(print) print_cores(core);

	delete[] peel_sequence;
	delete[] core;

	printf("\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
}

void Graph::core_decomposition_hindex(bool print) {
#ifdef NDEBUG
	printf("*** core_decomposition h-index (Release): %s ***\n", dir.c_str());
#else
	printf("*** core_decomposition h-index (Debug): %s ***\n", dir.c_str());
#endif

	read_graph_binary();

	Timer timer;
	ui *core = new ui[n]; //core[i] is the core value of vertex i

	core_decomposition_hindex(core);

	if(print) print_cores(core);

	delete[] core;

	printf("\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
}

void Graph::core_hierarchy(bool print) {
#ifdef NDEBUG
	printf("*** core_hierarchy (Release): %s ***\n", dir.c_str());
#else
	printf("*** core_hierarchy (Debug): %s ***\n", dir.c_str());
#endif

	read_graph_binary();

	Timer timer;
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n]; //core[i] is the core value of vertex i

	ui max_core = core_decomposition(peel_sequence, core);

#ifndef NDEBUG
	for(ui i = 1;i < n;i ++) assert(core[peel_sequence[i]] >= core[peel_sequence[i-1]]);
#endif

	UnionFind *uf = new UnionFind(n);
	uf->init(n);

	vector<pair<pair<ui,ui>, ui> > spt;

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]]) {
			ui v = edges[j];
			ui ru = uf->UF_find(u);
			ui rv = uf->UF_find(v);

			if(ru != rv) {
				spt.pb(make_pair(make_pair(u, v), core[u]));
				uf->UF_union(ru, rv);
			}
		}
		vis[u] = 1;
	}

	if(print) {
		FILE *f = Utility::open_file((dir + string("/core_spanning_tree.txt")).c_str(), "w");
		fprintf(f, "vertex vertex weight\n");
		for(ui i = spt.size();i > 0;i --) fprintf(f, "%u\t%u\t%u\n", spt[i-1].first.first, spt[i-1].first.second, spt[i-1].second);
		fclose(f);
	}

	delete[] vis;
	delete[] peel_sequence;
	delete[] core;

	printf("\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
}

// private member functions

void Graph::read_graph_binary() {
	printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
	FILE *f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(ui), 1, f);
	if(tt != sizeof(ui)) {
		printf("sizeof unsigned int is different: b_degree.bin(%u), machine(%lu)\n", tt, sizeof(ui));
		return ;
	}

	fread(&n, sizeof(ui), 1, f);
	fread(&m, sizeof(ui), 1, f);

	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *degree = new ui[n];
	fread(degree, sizeof(ui), n, f);

	fclose(f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	assert(sum == m);
#endif

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[n+1];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
		pstart[i+1] = pstart[i] + degree[i];
	}

	fclose(f);

	delete[] degree;

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

ui Graph::core_decomposition(ui *peel_sequence, ui *core) {
#ifndef NDEBUG
	printf("Start core decomposition\n");
#endif
	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1]-pstart[i];

	ui *rid = new ui[n];
	ui *id = peel_sequence;
	memset(id, 0, sizeof(ui)*n);
	for(ui i = 0;i < n;i ++) ++ id[degree[i]];
	for(ui i = 1;i < n;i ++) id[i] += id[i-1];

	for(ui i = 0;i < n;i ++) rid[i] = -- id[degree[i]];
	for(ui i = 0;i < n;i ++) id[rid[i]] = i;

	ui *degree_start = new ui[n+1];
	for(ui i = 0, j = 0;i <= n;i ++) {
		while(j < n&&degree[id[j]] < i) ++ j;
		degree_start[i] = j;
	}

	ui max_core = 0;
	for(ui i = 0;i < n;i ++) {
		ui u = id[i];
		assert(degree_start[degree[u]] == i);
		if(degree[u] > max_core) max_core = degree[u];
		core[u] = max_core;

		++ degree_start[degree[u]];
		if(degree[u] == 0) continue;

		degree_start[degree[u]-1] = degree_start[degree[u]];
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(rid[edges[j]] > i) {
			ui v = edges[j];
			ui pos1 = degree_start[degree[v]], pos2 = rid[v];
			swap(id[pos1], id[pos2]);
			rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
			++ degree_start[degree[v]];
			-- degree[v];
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) {
		ui cnt = 0;
		ui u = peel_sequence[i];
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(rid[edges[j]] > i) ++ cnt;
		assert(cnt == degree[u]);
	}
#endif

	delete[] degree;
	delete[] degree_start;
	delete[] rid;

#ifndef NDEBUG
	printf("Finished core decomposition\n");
#endif
	return max_core;
}

ui Graph::core_decomposition_hindex(ui *core) {
#ifndef NDEBUG
	printf("Start core decomposition\n");
#endif
	for(ui i = 0;i < n;i ++) core[i] = pstart[i+1]-pstart[i];

	ui *cnt = new ui[n];
	ui *q_1 = new ui[n];
	ui *q_2 = new ui[n];
	ui q_1_n = 0;
	for(ui i = 0;i < n;i ++) {
		cnt[i] = 0;
		for(ui j = pstart[i];j < pstart[i+1];j ++) if(core[edges[j]] >= core[i]) ++ cnt[i];
		if(cnt[i] < core[i]) q_1[q_1_n++] = i;
	}

	char *in_q = new char[n];
	memset(in_q, 0, sizeof(char)*n);

	ui *hindex = new ui[n];

#ifndef NDEBUG
	ui cc = 0;
	printf("iteration queue_size\n");
#endif
	while(q_1_n) {
#ifndef NDEBUG
		printf("%u %u\n", ++ cc, q_1_n);
#endif
		ui q_2_n = 0;
		for(ui i = 0;i < q_1_n;i ++) in_q[q_1[i]] = 1;
		for(ui i = 0;i < q_1_n;i ++) {
			ui u = q_1[i]; in_q[u] = 0;
			ui old = core[u];

			for(ui j = 0;j <= old;j ++) hindex[j] = 0;
			for(ui j = pstart[u];j < pstart[u+1];j ++) {
				ui val = core[edges[j]];
				if(val > old) ++ hindex[old];
				else ++ hindex[val];
			}
			for(ui j = old;j > 0;j --) {
				if(hindex[j] >= j) {
					core[u] = j;
					break;
				}
				hindex[j-1] += hindex[j];
			}

			//if(old != core[u]) {
				cnt[u] = 0;
				for(ui j = pstart[u];j < pstart[u+1];j ++) {
					ui v = edges[j];
					if(core[v] >= core[u]) ++ cnt[u];
					if(!in_q[v]&&core[u] < core[v]&&old >= core[v]) {
						-- cnt[v];
						if(cnt[v] +1 == core[v]) q_2[q_2_n++] = v;
					}
				}
			//}
		}

		q_1_n = q_2_n;
		swap(q_1, q_2);
	}

	ui max_core = 0;
	for(ui i = 0;i < n;i ++) if(core[i] > max_core) max_core = core[i];

	delete[] hindex;
	delete[] in_q;
	delete[] q_1;
	delete[] q_2;
	delete[] cnt;

#ifndef NDEBUG
	printf("Finished core decomposition\n");
#endif
	return max_core;
}

ui Graph::core_decomposition_linear_heap(ui alg, ui *peel_sequence, ui *core) {
#ifndef NDEBUG
	printf("Start core decomposition\n");
#endif
	if(alg != 1&&alg != 2) {
		printf("Parameter [alg] is wrong in core_decomposition_linear_heap\n");
		return 0;
	}

	ui *id_s = peel_sequence, *degree = core;
	for(ui i = 0;i < n;i ++) {
		id_s[i] = i;
		degree[i] = pstart[i+1]-pstart[i];
	}

	ui max_core = 0;
	if(alg == 1) {
		ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
		linear_heap->init(n, n-1, id_s, degree);
		memset(core, 0, sizeof(ui)*n);
		for(ui i = 0;i < n;i ++) {
			ui u, key;
			linear_heap->pop_min(u, key);
			if(key > max_core) max_core = key;
			peel_sequence[i] = u;
			core[u] = max_core;

			for(ui j = pstart[u];j < pstart[u+1];j ++) if(core[edges[j]] == 0) {
				linear_heap->decrement(edges[j]);
			}
		}
		delete linear_heap;
	}
	else {
		ArrayLinearHeap *linear_heap = new ArrayLinearHeap(n, n-1);
		linear_heap->init(n, n-1, id_s, degree);
		memset(core, 0, sizeof(ui)*n);
		for(ui i = 0;i < n;i ++) {
			ui u, key;
			linear_heap->pop_min(u, key);
			if(key > max_core) max_core = key;
			peel_sequence[i] = u;
			core[u] = max_core;

			for(ui j = pstart[u];j < pstart[u+1];j ++) if(core[edges[j]] == 0) {
				linear_heap->decrement(edges[j]);
			}
		}
		delete linear_heap;
	}

#ifndef NDEBUG
	printf("Finished core decomposition\n");
#endif
	return max_core;
}

void Graph::compute_core_gccs(const ui *seq, const ui *key) {
	for(ui i = 1;i < n;i ++) assert(key[seq[i]] >= key[seq[i-1]]);

	ui *cc_n = new ui[n];
	ui *cc_m = new ui[n];

	for(ui i = 0;i < n;i ++) {
		cc_n[i] = 1;
		cc_m[i] = 0;
	}

	FILE *f = Utility::open_file((dir + string("/core-gccs.txt")).c_str(), "w");
	fprintf(f, "core total_n total_m gcc_n gcc_m\n");
	ui t_n = 0, t_m = 0, max_n = 0, max_m = 0;

	UnionFind *uf = new UnionFind(n);
	uf->init();

	for(ui i = n;i > 0;i --) {
		ui u = seq[i-1];
		++ t_n;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(key[edges[j]] > key[u]||(key[edges[j]] == key[u]&&edges[j] > u)) {
			++ t_m;

			ui v = edges[j];
			ui tu = uf->UF_find(u);
			ui tv = uf->UF_find(v);

			ui r_u;
			if(tu == tv) {
				++ cc_m[tu];
				r_u = tu;
			}
			else {
				r_u = uf->UF_union(tu, tv);
				if(r_u == tu) {
					cc_n[tu] += cc_n[tv];
					cc_m[tu] += cc_m[tv] + 1;
				}
				else {
					assert(r_u == tv);
					cc_n[tv] += cc_n[tu];
					cc_m[tv] += cc_m[tu]+1;
				}
			}

			if(cc_n[r_u] > max_n) max_n = cc_n[r_u];
			if(cc_m[r_u] > max_m) max_m = cc_m[r_u];
		}
		if(i == 1||key[u] != key[seq[i-2]]) {
			ui low = 0;
			if(i > 1) low = key[seq[i-2]];
			for(ui j = key[u];j > low;j --) {
				fprintf(f, "%u %u %u %u %u\n", j, t_n, t_m, max_n, max_m);
			}
		}
	}

	fclose(f);

	delete uf;
	delete[] cc_n;
	delete[] cc_m;
}

void Graph::print_cores(const ui *core) {
	FILE *f = Utility::open_file((dir + string("/core_decomposition.txt")).c_str(), "w");
	fprintf(f, "vertex core_id\n");
	for(ui i = 0;i < n;i ++) fprintf(f, "%u\t%u\n", i, core[i]);
	fclose(f);
}
