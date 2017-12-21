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

/*
void Graph::parseDIMACS2binary() {
	vector<int> nodes;
	vector<pair<int,int> > edges;
	FILE *f = open_file((dir + string("/edges.txt")).c_str(), "r");

	char buf[1024];
	int a, b, c;
	for(ui i = 0;i < 8;i ++) {
		fgets(buf, 1024, f);
		printf("%s", buf);
	}
	while(fscanf(f, "%s%d%d", buf, &a, &b) ==3) {
		//printf("%d %d\n", a, b);
		if(a == b) continue;
		nodes.push_back(a);
		nodes.push_back(b);
		edges.push_back(make_pair(a,b));
		edges.push_back(make_pair(b,a));
	}

	fclose(f);

	write_graph_binary(nodes, edges);
}

void Graph::compress_byte() {
	string name = string(dir) + "/degree.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, fin);
	if(tt != sizeof(int)) {
		printf("the sizeof int is different edge_file(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, fin);

#ifdef _DEBUG_
	long long sum_degree = 0;
	for(ui i = 0;i < n;i ++) sum_degree += degree[i];
	if(sum_degree != m) {
		printf("WA #edges %lld m(%lld)\n", sum_degree, m);
		// return ;
	}
#endif

	fclose(fin);

	ui max_d = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > max_d) max_d = degree[i];
	printf("max_d: %u\n", max_d);

	name = string(dir) + "/adj.bin";
	fin = open_file(name.c_str(), "rb");

	ui *buf = new ui[max_d];

	unsigned char *obuf = new unsigned char[5*max_d];

	name = string(dir) + "/compress_b_adj.bin";
	FILE *fout = open_file(name.c_str(), "wb");

	name = string(dir) + "/compress_b_len.bin";
	FILE *flen = open_file(name.c_str(), "wb");

	int size = sizeof(int);
	fwrite(&size, sizeof(int), 1, flen);
	fwrite(&n, sizeof(int), 1, flen);
	fwrite(&m, sizeof(int), 1, flen);

	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(buf, sizeof(int), degree[i], fin);

		ui len = 0;

		if(degree[i] == 0) {
			fwrite(&len, sizeof(int), 1, flen);
			continue;
		}

		sort(buf, buf+degree[i]);

		encode_byte_compressed(buf, degree[i], obuf, len);

		fwrite(&len, sizeof(int), 1, flen);

		fwrite(obuf, sizeof(unsigned char), len, fout);
	}

	fclose(fin);
	fclose(fout);
	fclose(flen);

	delete[] degree;
	delete[] buf;
	delete[] obuf;
}

void Graph::uncompress_byte() {
	string name = string(dir) + "/b_compress_b_len.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int size;
	fread(&size, sizeof(int), 1, fin);
	if(size != sizeof(int)) {
		printf("WA sizeof(int)!\n");
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %d, m = %ld\n", n, m);

	ui *len = new ui[n];
	fread(len, sizeof(int), n, fin);

	fclose(fin);

	ui *degree = new ui[n];

	name = string(dir) + "/b_compress_b_adj.bin";
	fin = open_file(name.c_str(), "rb");

	unsigned char *buf = new unsigned char[5*n];
	ui *adj = new ui[n];

	name = string(dir) + "adj.bin";
	FILE *fout = open_file(name.c_str(), "wb");

	for(ui i = 0;i < n;i ++) {
		ui d = 0;

		if(len[i] > 0) {
			fread(buf, sizeof(unsigned char), len[i], fin);

			decode_byte_compressed(buf, len[i], adj, d);
		}

		degree[i] = d;

		if(d > 0) fwrite(adj, sizeof(int), d, fout);
	}

	fclose(fin);
	fclose(fout);

	name = string(dir) + "degree.bin";
	fout = open_file(name.c_str(), "wb");

	fwrite(&size, sizeof(int), 1, fout);
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m, sizeof(int), 1, fout);
	fwrite(degree, sizeof(int), n, fout);

	delete[] degree;
	delete[] len;
	delete[] buf;
	delete[] adj;
}

void Graph::compress_nibble() {
	string name = string(dir) + "/b_degree.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, fin);
	if(tt != sizeof(int)) {
		printf("the sizeof int is different edge_file(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, fin);

#ifdef _DEBUG_
	long long sum_degree = 0;
	for(ui i = 0;i < n;i ++) sum_degree += degree[i];
	if(sum_degree != m) {
		printf("WA #edges\n");
		return ;
	}
#endif

	fclose(fin);

	name = string(dir) + "/b_adj.bin";
	fin = open_file(name.c_str(), "rb");

	ui *buf = new ui[n];
	ui *len = new ui[n];

	unsigned char *obuf = new unsigned char[5*n];

	name = string(dir) + "/b_compress_n_adj.bin";
	FILE *fout = open_file(name.c_str(), "wb");

	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(buf, sizeof(int), degree[i], fin);

		len[i] = 0;

		if(degree[i] == 0) continue;

		encode_nibble_compressed(buf, degree[i], obuf, len[i]);

		fwrite(obuf, sizeof(unsigned char), len[i], fout);
	}

	fclose(fin);
	fclose(fout);

	name = string(dir) + "/b_compress_n_len.bin";
	fout = open_file(name.c_str(), "wb");

	int size = sizeof(int);
	fwrite(&size, sizeof(int), 1, fout);
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m, sizeof(int), 1, fout);
	fwrite(len, sizeof(int), n, fout);

	fclose(fout);

	delete[] degree;
	delete[] len;
	delete[] buf;
	delete[] obuf;
}

void Graph::uncompress_nibble() {
	string name = string(dir) + "/b_compress_n_len.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int size;
	fread(&size, sizeof(int), 1, fin);
	if(size != sizeof(int)) {
		printf("WA sizeof(int)!\n");
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *len = new ui[n];
	fread(len, sizeof(int), n, fin);

	fclose(fin);

	ui *degree = new ui[n];

	name = string(dir) + "adj.bin";
	FILE *fout = open_file(name.c_str(), "wb");

	name = string(dir) + "/b_compress_n_adj.bin";
	fin = open_file(name.c_str(), "rb");

	unsigned char *buf = new unsigned char[5*n];
	ui *adj = new ui[n];

	int cnt = 0;

	for(ui i = 0;i < n;i ++) {
		ui d = 0;

		if(len[i] > 0) {
			fread(buf, sizeof(unsigned char), len[i], fin);

			decode_nibble_compressed(buf, len[i], adj, d);
		}

		for(int j = 0;j < d;j ++) if(adj[j] == i) printf("self-loop %d!\n", cnt ++);

		int nd = 0;
		for(int j = 0;j < d;j ++) if(adj[j] != i) adj[nd ++] = adj[j];

		degree[i] = d;

		if(d > 0) fwrite(adj, sizeof(int), d, fout);
	}

	fclose(fin);
	fclose(fout);

	name = string(dir) + "degree.bin";
	fout = open_file(name.c_str(), "wb");

	ui tt = (ui)m;

	fwrite(&size, sizeof(int), 1, fout);
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&tt, sizeof(int), 1, fout);
	fwrite(degree, sizeof(int), n, fout);

	fclose(fout);

	delete[] degree;
	delete[] len;
	delete[] buf;
	delete[] adj;
}

void Graph::compress_and_orient(GraphStore igs, GraphStore ogs) {
	string name = string(dir) + "/b_degree.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int size;
	fread(&size, sizeof(int), 1, fin);
	if(size != sizeof(int)) {
		printf("WA sizeof(int)!\n");
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, fin);
	fclose(fin);

	ui *len = new ui[n];

	if(igs == byte_compressed) {
		name = string(dir) + "/b_compress_b_len.bin";
		fin = open_file(name.c_str(), "rb");
		fread(len, sizeof(int), 3, fin);
		fread(len, sizeof(int), n, fin);
		fclose(fin);

		name = string(dir) + "/b_compress_b_adj.bin";
		fin = open_file(name.c_str(), "rb");
	}
	else if(igs == nibble_compressed) {
		name = string(dir) + "/b_compress_n_len.bin";
		fin = open_file(name.c_str(), "rb");
		fread(len, sizeof(int), 3, fin);
		fread(len, sizeof(int), n, fin);
		fclose(fin);

		name = string(dir) + "/b_compress_n_adj.bin";
		fin = open_file(name.c_str(), "rb");
	}
	else {
		name = string(dir) + "/b_adj.bin";
		fin = open_file(name.c_str(), "rb");
	}

	FILE *fout = NULL;
	if(ogs == byte_compressed) {
		name = string(dir) + "/b_oriented_compress_b_adj.bin";
		fout = open_file(name.c_str(), "wb");
	}
	else if(ogs == nibble_compressed) {
		name = string(dir) + "/b_oriented_compress_n_adj.bin";
		fout = open_file(name.c_str(), "wb");
	}

	ui *adj = new ui[n];
	unsigned char *buf = new unsigned char[((long long)5)*n];
	for(ui i = 0;i < n;i ++) {
		ui tlen = 0;

		if(igs == uncompressed) fread(adj, sizeof(int), degree[i], fin);
		else if(igs == byte_compressed) {
			fread(buf, sizeof(unsigned char), len[i], fin);
			decode_byte_compressed(buf, len[i], adj, tlen);
#ifdef _DEBUG_
			if(tlen != degree[i]) printf("WA in decode\n");
#endif
		}
		else if(igs == nibble_compressed) {
			fread(buf, sizeof(unsigned char), len[i], fin);
			decode_nibble_compressed(buf, len[i], adj, tlen);
#ifdef _DEBUG_
			if(tlen != degree[i]) printf("WA in decode\n");
#endif
		}
		tlen = 0;
		for(ui j = 0;j < degree[i];j ++) {
			if(degree[adj[j]] > degree[i]||(degree[adj[j]] == degree[i]&&adj[j] > i)) {
				adj[tlen ++] = adj[j];
			}
		}

		if(ogs == byte_compressed) encode_byte_compressed(adj, tlen, buf, len[i]);
		else if(ogs == nibble_compressed) encode_nibble_compressed(adj, tlen, buf, len[i]);

		fwrite(buf, sizeof(unsigned char), len[i], fout);
	}

	fclose(fin);
	fclose(fout);

	if(ogs == byte_compressed) {
		name = string(dir) + "/b_oriented_compress_b_len.bin";
		fout = open_file(name.c_str(), "wb");
	}
	else if(ogs == nibble_compressed) {
		name = string(dir) + "/b_oriented_compress_n_len.bin";
		fout = open_file(name.c_str(), "wb");
	}

	m /= 2;

	fwrite(&size, sizeof(int), 1, fout);
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m, sizeof(int), 1, fout);
	fwrite(len, sizeof(int), n, fout);
	fclose(fout);

	delete[] buf;
	delete[] adj;
	delete[] degree;
	delete[] len;
}

void Graph::decode_byte_compressed(unsigned char *buf, ui len, ui *adj, ui &tlen) {
	if(len < 4) {
		tlen = 0;
		return ;
	}

	adj[0] = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];

	ui j = 4;
	tlen = 1;

	while(j < len) {
		ui dif = 0;
		while(buf[j]&(1<<7)) {
			dif += (buf[j]^(1<<7));
			dif <<= 7;
			++ j;
		}

		dif += buf[j];
		++ j;

		adj[tlen] = adj[tlen-1] + dif;
		++ tlen;
	}
}

void Graph::decode_nibble_compressed(unsigned char *buf, ui len, ui *adj, ui &tlen) {
	if(len < 4) {
		tlen = 0;
		return ;
	}

	adj[0] = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];

	ui j = 4;
	tlen = 1;
	int b = 8;

	while(j < len) {
		ui dif = 0;
		while(true){
			if(b == 8) {
				dif += ((buf[j]>>4)&((1<<3)-1));
#ifdef _DEBUG_
				if(((buf[j]>>4)&((1<<3)-1)) < 0) printf("WA1\n");
#endif
				b -= 4;

				if((buf[j]&(1<<7)) == 0) break;

				dif <<= 3;
			}
			else if(b == 4) {
				dif += (buf[j]&((1<<3)-1));
				++ j; b = 8;

				if((buf[j-1]&(1<<3)) == 0) break;

				dif <<= 3;
			}
		}

		if(dif == 0) break;

		adj[tlen] = adj[tlen-1] + dif;
		++ tlen;
	}
}

void Graph::encode_byte_compressed(ui *adj, ui tlen, unsigned char *obuf, ui &len) {
	if(tlen == 0) {
		len = 0;
		return ;
	}

	ui id = adj[0];
	obuf[0] = id>>24; id &= (1<<24)-1;
	obuf[1] = id>>16; id &= (1<<16)-1;
	obuf[2] = id>>8; id &= (1<<8)-1;
	obuf[3] = id;

	long long xlen = 4;

	for(ui j = 1;j < tlen;j ++) {
		ui dif = adj[j] - adj[j-1];
		if(dif <= 0) {
			printf("WA! not sorted!\n");
			return ;
		}

		int k = 28;
		while(k >= 0&&((dif>>k)==0)) k -= 7;

		while(k >= 0) {
			obuf[xlen] = dif>>k;
			if(k > 0) obuf[xlen] |= (1<<7);
			dif &= (1<<k)-1;

			++ xlen;
			k -= 7;
		}
	}
	if(xlen >= (((long long)1)<<32)) printf("length of an adjacency list too long\n");
	len = xlen;
}

void Graph::encode_nibble_compressed(ui *adj, ui tlen, unsigned char *obuf, ui &len) {
	if(tlen == 0) {
		len = 0;
		return ;
	}
	ui id = adj[0];
	obuf[0] = id>>24; id &= (1<<24)-1;
	obuf[1] = id>>16; id &= (1<<16)-1;
	obuf[2] = id>>8; id &= (1<<8)-1;
	obuf[3] = id;

	long long xlen = 3;
	int b = 0;

	for(ui j = 1;j < tlen;j ++) {
		ui dif = adj[j] - adj[j-1];
		if(dif <= 0) {
			printf("WA! not sorted!\n");
			return ;
		}

		int k = 30;
		while(k >= 0&&((dif>>k)==0)) k -= 3;

		while(k >= 0) {
			if((dif>>k) >= 8) printf("WA encode\n");

			if(b == 0) {
				++ xlen;

				obuf[xlen] = (dif>>k)<<4;
				if(k > 0) obuf[xlen] |= (1<<7);
				b = 4;
			}
			else {
				obuf[xlen] |= (dif>>k);
				if(k > 0) obuf[xlen] |= (1<<3);
				b = 0;
			}

			dif &= (1<<k)-1;

			k -= 3;
		}
	}

	++ xlen;
	if(xlen >= (((long long)1)<<32)) printf("length of an adjacency list too long\n");
	len = xlen;
}

void Graph::trilist_byte() {
#ifdef _LINUX_
	struct timeval start;
	gettimeofday(&start, NULL);
#else
	int start = clock();
#endif

	string name = dir + "/b_oriented_compress_b_len.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, fin);
	if(tt != sizeof(int)) {
		printf("the sizeof int is different edge_file(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *ps = new ui[n+1];
	fread(ps+1, sizeof(int), n, fin);
	fclose(fin);

	ps[0] = 0;
	for(ui i = 0;i < n;i ++) {
#ifdef _DEBUG_
		if(((long long)ps[i]) + ps[i+1] >= (((long long)1)<<32)) {
			printf("Too long\n");
			return ;
		}
#endif
		ps[i+1] += ps[i];
	}

	unsigned char *codes = new unsigned char[ps[n]];

	name = dir + "/b_oriented_compress_b_adj.bin";
	fin = open_file(name.c_str(), "rb");

	fread(codes, sizeof(unsigned char), ps[n], fin);
	fclose(fin);

#ifdef _LINUX_
	struct timeval end1;
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#else
	int end1 = clock();
#endif

	char *adj = new char[n>>3];
	memset(adj, 0, sizeof(char)*(n>>3));

	long long res = 0;

	for(ui i = 0;i < n;i ++) {
		unsigned char *buf = codes + ps[i];
		ui len = ps[i+1] - ps[i];
		if(len == 0) continue;

		ui id = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];
		adj[id>>3] |= (1<<(id&7));
		for(ui j = 4;j < len;) {
			ui dif = 0;
			while(buf[j]&(1<<7)) {
				dif += (buf[j]^(1<<7));
				dif <<= 7;
				++ j;
			}

			dif += buf[j];
			++ j;

			id += dif;
			adj[id>>3] |= (1<<(id&7));
		}

		id = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];
		ui j = 4;
		while(true) {
			unsigned char *tbuf = codes + ps[id];
			ui tlen = ps[id+1] - ps[id];

			ui tid;
			if(tlen > 0) {
				tid = (((((((ui)tbuf[0])<<8)+tbuf[1])<<8)+tbuf[2])<<8) + tbuf[3];
				if(adj[tid>>3]&(1<<(tid&7))) ++ res;
			}

			for(ui k = 4;k < tlen;) {
				ui dif = 0;
				while(tbuf[k]&(1<<7)) {
					dif += (tbuf[k]^(1<<7));
					dif <<= 7;
					++ k;
				}

				dif += tbuf[k];
				++ k;

				tid += dif;
				if(adj[tid>>3]&(1<<(tid&7))) ++ res;
			}

			if(j >= len) break;

			ui dif = 0;
			while(buf[j]&(1<<7)) {
				dif += (buf[j]^(1<<7));
				dif <<= 7;
				++ j;
			}

			dif += buf[j];
			++ j;

			id += dif;
		}

		id = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];
		adj[id>>3] = 0;
		for(ui j = 4;j < len;) {
			ui dif = 0;
			while(buf[j]&(1<<7)) {
				dif += (buf[j]^(1<<7));
				dif <<= 7;
				++ j;
			}

			dif += buf[j];
			++ j;

			id += dif;
			adj[id>>3] = 0;
		}
	}

	printf("#Triangles by trilist-byte: %lld\n", res);

#ifdef _LINUX_
	struct timeval end;
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - end1.tv_sec;
	useconds = end.tv_usec - end1.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("IO time: %lld\nTriangle time: %lld\n", mtime1, mtime);
#else
	int end = clock();

	printf("IO time: %d\nTriangle time: %d\n", end1-start,end-end1);
#endif

	delete[] ps;
	delete[] codes;
	delete[] adj;
}

void Graph::trilist_byte_plain() {
#ifdef _LINUX_
	struct timeval start;
	gettimeofday(&start, NULL);
#else
	int start = clock();
#endif

	string name = dir + "/b_oriented_compress_b_len.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, fin);
	if(tt != sizeof(int)) {
		printf("the sizeof int is different edge_file(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *ps = new ui[n+1];
	ui *len = new ui[n];
	ui *edges = new ui[m];

	fread(len, sizeof(int), n, fin);
	fclose(fin);

	unsigned char *codes = new unsigned char[5*n];

	ui *degree = new ui[n];
	name = dir + "/b_degree.bin";
	fin = open_file(name.c_str(), "rb");
	fread(degree, sizeof(int), 3, fin);
	fread(degree, sizeof(int), n, fin);
	fclose(fin);

	name = dir + "/b_adj.bin";
	FILE *fedge = open_file(name.c_str(), "rb");

	name = dir + "/b_oriented_compress_b_adj.bin";
	fin = open_file(name.c_str(), "rb");

	ui *buf = new ui[n];
	ps[0] = 0;
	for(ui i = 0;i < n;i ++) {
		fread(codes, sizeof(unsigned char), len[i], fin);
		decode_byte_compressed(codes, len[i], edges+ps[i], ps[i+1]);
		ps[i+1] += ps[i];

		fread(buf, sizeof(int), degree[i], fedge);
		ui tlen = 0;
		for(ui j = 0;j < degree[i];j ++) if(degree[buf[j]] > degree[i]||(degree[buf[j]] == degree[i]&&buf[j] > i)) {
			buf[tlen ++] = buf[j];
		}

		if(tlen != ps[i+1] - ps[i]) printf("Len WA %u %u\n", tlen, ps[i+1]-ps[i]);
		for(ui j = 0;j < tlen;j ++) if(buf[j] != edges[ps[i]+j]) printf("ID WA\n");
	}
	fclose(fin);
	fclose(fedge);

	delete[] buf;

#ifdef _LINUX_
	struct timeval end1;
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#else
	int end1 = clock();
#endif
	char *adj = new char[n];
	memset(adj, 0, sizeof(char)*n);

	long long res = 0;

	for(ui i = 0;i < n;i ++) {
		for(ui j = ps[i];j < ps[i+1];j ++) adj[edges[j]] = 1;

		for(ui j = ps[i];j < ps[i+1];j ++) {
			int v = edges[j];
			for(ui k = ps[v];k < ps[v+1];k ++) {
				if(adj[edges[k]]) ++ res;
			}
		}

		for(ui j = ps[i];j < ps[i+1];j ++) adj[edges[j]] = 0;
	}

	printf("#Triangles by trilist-byte: %lld\n", res);

#ifdef _LINUX_
	struct timeval end;
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - end1.tv_sec;
	useconds = end.tv_usec - end1.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("IO time: %lld\nTriangle time: %lld\n", mtime1, mtime);
#else
	int end = clock();

	printf("IO time: %d\nTriangle time: %d\n", end1-start,end-end1);
#endif

	delete[] degree;
	delete[] len;
	delete[] ps;
	delete[] codes;
	delete[] adj;
}

void Graph::trilist_nibble() {
#ifdef _LINUX_
	struct timeval start;
	gettimeofday(&start, NULL);
#else
	int start = clock();
#endif

	string name = dir + "/b_oriented_compress_n_len.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, fin);
	if(tt != sizeof(int)) {
		printf("the sizeof int is different edge_file(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *ps = new ui[n+1];
	fread(ps+1, sizeof(int), n, fin);
	fclose(fin);

	ps[0] = 0;
	for(ui i = 0;i < n;i ++) {
#ifdef _DEBUG_
		if(((long long)ps[i]) + ps[i+1] >= (((long long)1)<<32)) {
			printf("Too long\n");
			return ;
		}
#endif
		ps[i+1] += ps[i];
	}

	unsigned char *codes = new unsigned char[ps[n]];

	name = dir + "/b_oriented_compress_n_adj.bin";
	fin = open_file(name.c_str(), "rb");

	fread(codes, sizeof(unsigned char), ps[n], fin);
	fclose(fin);

#ifdef _LINUX_
	struct timeval end1;
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#else
	int end1 = clock();
#endif

	char *adj = new char[n>>3];
	memset(adj, 0, sizeof(char)*(n>>3));

	long long res = 0;

	for(ui i = 0;i < n;i ++) {
		unsigned char *buf = codes + ps[i];
		ui len = ps[i+1] - ps[i];
		if(len == 0) continue;

		// set arrayset
		ui id = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];
		adj[id>>3] |= (1<<(id&7));

		for(ui j = 4, b = 8;j < len;) {
			ui dif = 0;
			while(true){
				if(b == 8) {
					dif += ((buf[j]>>4)&7);
					b -= 4;
					if((buf[j]&(1<<7)) == 0) break;
					dif <<= 3;
				}
				else if(b == 4) {
					dif += (buf[j]&7);
					++ j; b = 8;
					if((buf[j-1]&8) == 0) break;
					dif <<= 3;
				}
			}

			if(dif == 0) break;

			id += dif;
			adj[id>>3] |= (1<<(id&7));
		}

		//do intersection
		id = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];
		ui j = 4, b= 8;
		while(true) {
			unsigned char *tbuf = codes + ps[id];
			ui tlen = ps[id+1] - ps[id];

			ui tid;
			if(tlen > 0) {
				tid = (((((((ui)tbuf[0])<<8)+tbuf[1])<<8)+tbuf[2])<<8) + tbuf[3];
				if(adj[tid>>3]&(1<<(tid&7))) ++ res;
			}

			for(ui k = 4, bb = 8;k < tlen;) {
				ui dif = 0;
				while(true){
					if(bb == 8) {
						dif += ((tbuf[k]>>4)&7);
						bb -= 4;
						if((tbuf[k]&(1<<7)) == 0) break;
						dif <<= 3;
					}
					else if(bb == 4) {
						dif += (tbuf[k]&7);
						++ k; bb = 8;
						if((tbuf[k-1]&8) == 0) break;
						dif <<= 3;
					}
				}

				if(dif == 0) break;

				tid += dif;
				if(adj[tid>>3]&(1<<(tid&7))) ++ res;
			}

			if(j >= len) break;

			ui dif = 0;
			while(true){
				if(b == 8) {
					dif += ((buf[j]>>4)&7);
					b -= 4;
					if((buf[j]&(1<<7)) == 0) break;
					dif <<= 3;
				}
				else if(b == 4) {
					dif += (buf[j]&7);
					++ j; b = 8;
					if((buf[j-1]&8) == 0) break;
					dif <<= 3;
				}
			}

			if(dif == 0) break;

			id += dif;
		}

		//  reset to 0
		id = (((((((ui)buf[0])<<8)+buf[1])<<8)+buf[2])<<8) + buf[3];
		adj[id>>3] = 0;

		for(ui j = 4, b = 8;j < len;) {
			ui dif = 0;
			while(true){
				if(b == 8) {
					dif += ((buf[j]>>4)&7);
					b -= 4;
					if((buf[j]&(1<<7)) == 0) break;
					dif <<= 3;
				}
				else if(b == 4) {
					dif += (buf[j]&7);
					++ j; b = 8;
					if((buf[j-1]&8) == 0) break;
					dif <<= 3;
				}
			}

			if(dif == 0) break;

			id += dif;
			adj[id>>3] = 0;
		}
	}

	printf("#Triangles by trilist-nibble: %lld\n", res);

#ifdef _LINUX_
	struct timeval end;
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - end1.tv_sec;
	useconds = end.tv_usec - end1.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("IO time: %lld\nTriangle time: %lld\n", mtime1, mtime);
#else
	int end = clock();

	printf("IO time: %d\nTriangle time: %d\n", end1-start,end-end1);
#endif

	delete[] ps;
	delete[] codes;
	delete[] adj;
}

void Graph::trilist_chiba() {
#ifdef _LINUX_
	struct timeval start;
	gettimeofday(&start, NULL);
#else
	int start = clock();
#endif

	string name = dir + "/b_degree.bin";
	FILE *fin = open_file(name.c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, fin);
	if(tt != sizeof(int)) {
		printf("the sizeof int is different edge_file(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	printf("n = %u, m = %ld\n", n, m);

	ui *pstart = new ui[n+1];
	fread(pstart+1, sizeof(int), n, fin);
	fclose(fin);

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) pstart[i+1] += pstart[i];

	ui *edges = new ui[m];
	name = dir + "/b_adj.bin";
	fin = open_file(name.c_str(), "rb");
	for(ui i = 0;i < n;i ++) {
		fread(edges+pstart[i], sizeof(int), pstart[i+1]-pstart[i], fin);
	}
	fclose(fin);

#ifdef _LINUX_
	struct timeval end1;
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#else
	int end1 = clock();
#endif

	ui *bin_head = new ui[n];
	ui *bin_next = new ui[n];
	for(ui i = 0;i < n;i ++) bin_head[i] = n;

	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		bin_next[i] = bin_head[d];
		bin_head[d] = i;
	}

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	char *adj = new char[n];
	memset(adj, 0, sizeof(char)*n);

	long long res = 0;

	ui *vv = new ui[n];
	ui vv_n;

	int max_d = n-1;
	for(ui i = 0;i < n;i ++) {
		while(bin_head[max_d] == n) -- max_d;

		ui u = bin_head[max_d];
		bin_head[max_d] = bin_next[u];

		vv_n = 0;

		for(ui j = pstart[u];j < pstart[u+1];j ++) {
			if(!vis[edges[j]]) {
				vv[vv_n ++] = edges[j];
				adj[edges[j]] = 1;
			}
		}

		for(ui j = 0;j < vv_n;j ++) {
			int v = vv[j];

			for(ui k = pstart[v];k < pstart[v+1]&&edges[k] < v;k ++) {
				if(adj[edges[k]]) ++ res;
			}
		}

		for(ui j = 0;j < vv_n;j ++) adj[vv[j]] = 0;

		vis[u] = 1;
	}

	printf("#Triangles by chiba: %lld\n", res);

#ifdef _LINUX_
	struct timeval end;
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - end1.tv_sec;
	useconds = end.tv_usec - end1.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("IO time: %lld\nTriangle time: %lld\n", mtime1, mtime);
#else
	int end = clock();

	printf("IO time: %d\nTriangle time: %d\n", end1-start,end-end1);
#endif

	delete[] pstart;
	delete[] edges;
	delete[] bin_head;
	delete[] bin_next;
	delete[] vis;
	delete[] vv;
	delete[] adj;
}

void Graph::triangle_count_for_edges() {
	ui *deg = new ui[n];
		for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

		ui *adj = new ui[n];
		memset(adj, 0, sizeof(ui)*n);

		ui *similar = new ui[m];
		memset(similar, 0, sizeof(ui)*m);

		ui *pend = new ui[n];

		for(ui i = 0;i < n;i ++) {
			ui &end = pend[i] = pstart[i];
			ui j = pstart[i+1];
			while(true) {
				while(end < j&&(deg[edges[end]] < deg[i]||(deg[edges[end]]==deg[i]&&edges[end]<i))) ++ end;
				while(j > end&&(deg[edges[j-1]] > deg[i]||(deg[edges[j-1]]==deg[i]&&edges[j-1]>i))) -- j;
				if(end >= j) break;
				swap(edges[end], edges[j-1]);
			}
			sort(edges+pend[i], edges+pstart[i+1]);
		}

		for(ui u = 0;u < n;u ++) {
			for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

			for(ui j = pstart[u];j < pend[u];j ++) {
				ui v = edges[j];

				for(ui k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) {
					++ similar[j];
					++ similar[k];
					++ similar[adj[edges[k]] - 1];
				}
			}

			for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
		}
}*/
