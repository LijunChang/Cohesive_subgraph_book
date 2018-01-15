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

void Graph::k_edge_connected_component(ui K, bool print) {
#ifdef NDEBUG
	printf("*** k_edge_connected_component (Release): %s, %u ***\n", dir.c_str(), K);
#else
	printf("*** k_edge_connected_component (Debug): %s, %u ***\n", dir.c_str(), K);
#endif

	if(K < 2) {
		printf("K must be at least 2!\n");
		return ;
	}

	read_graph_binary();

	Timer timer;

	ui *pend = new ui[n];
	for(ui i = 0;i < n;i ++) pend[i] = pstart[i+1];

	char *computed = new char[n];
	memset(computed, 0, sizeof(char)*n);

	ui *degree = new ui[n];
	ui *Q = new ui[n];
	ui Q_n = 0;

	// k-core-based pruning
	for(ui i = 0;i < n;i ++) {
		degree[i] = pend[i] - pstart[i];
		if(degree[i] < K) {
			Q[Q_n ++] = i;
			computed[i] = 1;
		}
	}
	k_core_prune(K, Q, Q_n, computed, degree, pend);

	ui *ids = new ui[n];
	ui *cstart = new ui[n+1];
	ui c_n = 0;
	cstart[0] = 0;

	Edge **graph_head = new Edge*[n];
	Edge *graph_edges = new Edge[m];
	for(ui i = 0;i < n;i ++) graph_head[i] = nullptr;

	ui *sv_next = new ui[n];
	ui *sv_last = new ui[n];

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	ListLinearHeap *heap = new ListLinearHeap(n, K);

	ui *keys = new ui[n];

	for(ui i = 0;i < n;) {
		if(computed[i]) {
			++ i;
			continue;
		}

#ifndef NDEBUG
		for(ui j = 0;j < n;j ++) assert(vis[j] == 0&&graph_head[j] == nullptr);
#endif
		construct_pgraph(i, Q, vis, computed, pend, graph_head, graph_edges, sv_next, sv_last) ;
		ui new_cn = decomposition(i, K, cstart, ids, c_n, Q, keys, vis, graph_head, graph_edges, heap, sv_next, sv_last);
		//printf("c_n: %u, new_cn: %u\n", c_n, new_cn);

		/*for(ui j = c_n;j < new_cn;j ++) {
			printf("ids of %u:", j);
			for(ui k = cstart[j];k < cstart[j+1];k ++) printf(" %u", ids[k]);
			printf("\n");
		}*/

		if(new_cn == c_n + 1) {
			for(ui j = cstart[c_n];j < cstart[c_n+1];j ++) computed[ids[j]] = 1;
			++ c_n;
		}
		else remove_inter_edges(K, c_n, new_cn, cstart, ids, keys, pend, Q, computed, degree); // use array keys for component id
	}

	if(print) print_kecc(K, c_n, cstart, ids);

	delete[] keys;
	delete heap;
	delete[] vis;
	delete[] sv_next;
	delete[] sv_last;
	delete[] graph_head;
	delete[] graph_edges;
	delete[] cstart;
	delete[] ids;
	delete[] Q;
	delete[] degree;
	delete[] computed;
	delete[] pend;

	printf("\tTotal processing time excluding I/O: %s (microseconds)\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
}

void Graph::k_edge_connected_component_space(ui K, bool print) {
#ifdef NDEBUG
	printf("*** k_edge_connected_component_space (Release): %s, %u ***\n", dir.c_str(), K);
#else
	printf("*** k_edge_connected_component_space (Debug): %s, %u ***\n", dir.c_str(), K);
#endif

	if(K < 2) {
		printf("K must be at least 2!\n");
		return ;
	}

	read_graph_binary();

	Timer timer;

	ui *pend = new ui[n];
	for(ui i = 0;i < n;i ++) pend[i] = pstart[i+1];

	char *computed = new char[n];
	ui *Q = new ui[n], Q_n = 0;
	ui *degree = new ui[n];
	// k-core-based pruning
	memset(computed, 0, sizeof(char)*n);
	for(ui i = 0;i < n;i ++) {
		degree[i] = pend[i] - pstart[i];
		if(degree[i] < K) {
			Q[Q_n ++] = i;
			computed[i] = 1;
		}
	}
	k_core_prune(K, Q, Q_n, computed, degree, pend);

	UnionFind *UF = new UnionFind(n);
	ui *representative = new ui[n];
	ui *sv_next = new ui[n];
	ui *sv_last = new ui[n];
	ui *adj_next = new ui[n];
	ui *adj_last = new ui[n];

	ui *ids = new ui[n];
	ui *cstart = new ui[n+1];
	ui c_n = 0;
	cstart[0] = 0;

	ListLinearHeap *heap = new ListLinearHeap(n, K);
	ui *keys = new ui[n];

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	ui *pend_local = new ui[n];

	for(ui i = 0;i < n;) {
		if(computed[i]) {
			++ i;
			continue;
		}

		initialize_pgraph(i, Q, vis, computed, pend, pend_local, UF, representative, sv_next, sv_last, adj_next, adj_last);
		ui new_cn = decomposition(i, K, cstart, ids, c_n, Q, keys, vis, pend_local, heap, sv_next, sv_last, adj_next, adj_last, UF, representative);
		if(new_cn == c_n + 1) {
			for(ui j = cstart[c_n];j < cstart[c_n+1];j ++) computed[ids[j]] = 1;
			++ c_n;
		}
		else remove_inter_edges(K, c_n, new_cn, cstart, ids, pend, Q, computed, degree, UF);
	}

	if(print) print_kecc(K, c_n, cstart, ids);

	delete[] pend_local;
	delete[] vis;
	delete[] keys;
	delete[] pend;

	delete[] ids;
	delete[] cstart;

	delete[] representative;
	delete UF;
	delete[] sv_next;
	delete[] sv_last;
	delete[] adj_next;
	delete[] adj_last;

	delete[] computed;
	delete[] degree;
	delete[] Q;

	printf("\tTotal processing time excluding I/O: %s (microseconds)\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
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

void Graph::k_core_prune(ui K, ui *Q, ui Q_n, char *computed, ui *degree, ui *pend) {
	for(ui i = 0;i < Q_n;i ++) {
		ui u = Q[i];
		for(ui j = pstart[u];j < pend[u];j ++) if(!computed[edges[j]]) {
			ui v = edges[j];
			if(degree[v] == K) {
				Q[Q_n ++] = v;
				computed[v] = 1;
			}
			-- degree[v];
		}
	}
#ifndef NDEBUG
	printf("k_core based pruning removed %u vertices!\n", Q_n);
#endif
}

void Graph::construct_pgraph(ui s, ui *Q, char *vis, char *computed, ui *pend, Edge **graph_head, Edge *graph_edges, ui *sv_next, ui *sv_last) {
	ui cnt = 0;
	Q[0] = s; vis[s] = 1;
	ui Q_n = 1;
	assert(!computed[s]);

	for(ui i = 0;i < Q_n;i ++) {
		ui u = Q[i];
		sv_next[u] = sv_last[u] = u;
		for(ui j = pstart[u];j < pend[u];) {
			ui v = edges[j];
			if(computed[v]) swap(edges[j], edges[-- pend[u]]);
			else {
				if(!vis[v]) {
					Q[Q_n ++] = v;
					vis[v] = 1;
				}

				if(v > u) {
					graph_edges[cnt].vertex = v;
					graph_edges[cnt].reverse = &graph_edges[cnt+1];
					add_edge(u, &graph_edges[cnt], graph_head);
					++ cnt;

					graph_edges[cnt].vertex = u;
					graph_edges[cnt].reverse = &graph_edges[cnt-1];
					add_edge(v, &graph_edges[cnt], graph_head);
					++ cnt;
				}
				++ j;
			}
		}
	}
	for(ui i = 0;i < Q_n;i ++) vis[Q[i]] = 0;

	/*for(ui i = 0;i < Q_n;i ++) {
		printf("neighbors of %u:", Q[i]);
		for(Edge *e = graph_head[Q[i]];e != nullptr;e = e->next) printf(" %u", e->vertex);
		printf("\n");
	}*/
}

void Graph::initialize_pgraph(ui s,ui *Q, char *vis, char *computed, ui *pend, ui *pend_local, UnionFind *UF, ui *representative, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last) {
	ui Q_n = 1; Q[0] = s; vis[s] = 1;
	for(ui j = 0;j < Q_n;j ++) {
		ui u = Q[j];
		for(ui k = pstart[u];k < pend[u];) {
			ui v = edges[k];
			if(computed[v]) swap(edges[k], edges[-- pend[u]]);
			else {
				if(!vis[v]) {
					Q[Q_n ++] = v;
					vis[v] = 1;
				}
				++ k;
			}
		}
	}
	UF->init(Q, Q_n);
	for(ui j = 0;j < Q_n;j ++) {
		ui u = Q[j];
		vis[u] = 0;
		representative[u] = u;
		sv_next[u] = sv_last[u] = adj_next[u] = adj_last[u] = u;
		pend_local[u] = pend[u];
	}
}

ui Graph::decomposition(ui s, ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, Edge **graph_head, Edge *graph_edges, ListLinearHeap *heap, ui *sv_next, ui *sv_last) {
	while(graph_head[s] != nullptr) {
		heap->init(0, K, nullptr, nullptr);
		heap->insert(s, 0);
		ui u, key, Q_n = 0;
		while(heap->pop_max(u, key)) {
			Q[Q_n ++] = u; vis[u] = 1; //vis[u] = 1 means u is in Q, and vis[u] = 2 means u is in heap
										//vis[u] = 3 means u is in Q between Q_n and new_Qn
			keys[u] = key;

			ui new_Qn = Q_n;
			for(ui j = Q_n-1;j < new_Qn;j ++) {
				ui v = Q[j];

				for(Edge *e = graph_head[v];e != nullptr;e = e->next) if(vis[e->vertex] != 1) {
					ui w = e->vertex;

					if(vis[w] == 3) {
						++ keys[w];
						continue;
					}

					if(vis[w] == 2) key = heap->remove(w);
					else key = 0;
					assert(key < K);

					++ key;
					if(key >= K) {
						Q[new_Qn ++] = w;
						keys[w] = key;
						vis[w] = 3;
					}
					else {
						heap->insert(w, key);
						vis[w] = 2;
					}
				}

				if(v == u) continue;

				// contract u and v
				vis[v] = 0;
				keys[u] += keys[v];
				merge(graph_head, graph_edges, u, v, keys, sv_next, sv_last);
			}
		}

		/*printf("before\n");
		for(ui i = 0;i < Q_n;i ++) {
			printf("neigbors of %u:", Q[i]);
			for(Edge *e = graph_head[Q[i]];e != nullptr;e = e->next) printf(" %u", e->vertex);
			printf("\n");
		}*/

		while(Q_n > 0&&keys[Q[Q_n-1]] < K) {
			ui u = Q[-- Q_n];
			vis[u] = 0;
			//printf("after remove %u\n", u);

			for(Edge *e = graph_head[u];e != nullptr;e = e->next) delete_edge(e->vertex, e->reverse, graph_head);
			graph_head[u] = nullptr;

			ui pos = cstart[c_n];
			ids[pos ++] = u;
			while(sv_next[u] != u) {
				u = sv_next[u];
				ids[pos ++] = u;
			}
			cstart[++ c_n] = pos;

			/*for(ui i = 0;i < Q_n;i ++) {
				printf("neigbors of %u:", Q[i]);
				for(Edge *e = graph_head[Q[i]];e != nullptr;e = e->next) printf(" %u", e->vertex);
				printf("\n");
			}*/
		}

		for(ui i = 0;i < Q_n;i ++) vis[Q[i]] = 0;
	}
	return c_n;
}

ui Graph::decomposition(const ui s, const ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, ui *pend_local, ListLinearHeap *heap, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last, UnionFind *UF, ui *representative) {
	ui old_cn = c_n;
	while(true) {
		heap->init(0, K, nullptr, nullptr);
		heap->insert(s, 0);
		ui u, key, Q_n = 0;
		//printf("start an iteration\n");
		while(heap->pop_max(u, key)) {
			Q[Q_n ++] = u; vis[u] = 1; //vis[u] = 1 means u is in Q, and vis[u] = 2 means u is in heap
										//vis[u] = 3 means u is in Q between Q_n and new_Qn
										//vis[u] = 4 means u is removed from the pgraph
			keys[u] = key;
			//printf(" [%u,%u]\n", u, key);

			ui new_Qn = Q_n;
			for(ui i = Q_n-1;i < new_Qn;i ++) {
				ui v = Q[i];
				ui pre = v, tv = v;

				while(true) {
					assert(representative[UF->UF_find(v)] == v);
					for(ui j = pstart[tv];j < pend_local[tv];) {
						ui w = representative[UF->UF_find(edges[j])];
						if(vis[w] == 4||w == v) {
							swap(edges[j], edges[-- pend_local[tv]]);
							continue;
						}
						++ j;
						if(vis[w] != 1) {
							if(vis[w] == 3) {
								++ keys[w];
								continue;
							}

							if(vis[w] == 2) key = heap->remove(w);
							else key = 0;
							assert(key < K);

							++ key;
							if(key >= K) {
								Q[new_Qn ++] = w;
								keys[w] = key;
								vis[w] = 3;
							}
							else {
								heap->insert(w, key);
								vis[w] = 2;
							}
						}
					}
					ui next = adj_next[tv];
					if(pend_local[tv] > pstart[tv]) {
						adj_next[pre] = tv;
						pre = tv;
					}
					if(next == tv) break;
					else tv = next;
				}
				adj_last[v] = pre;
				adj_next[pre] = pre;

				if(v == u) continue;

				// contract u and v
				keys[u] += keys[v];
				//printf("merged %u and %u\n", u, v);
				tv = v;
				while(true) {
					for(ui j = pstart[tv];j < pend_local[tv];j ++) if(representative[UF->UF_find(edges[j])] == u) -- keys[u];
					if(adj_next[tv] == tv) break;
					else tv = adj_next[tv];
				}
				sv_next[sv_last[u]] = v;
				sv_last[u] = sv_last[v];
				adj_next[adj_last[u]] = v;
				adj_last[u] = adj_last[v];
				representative[UF->UF_union(u,v)] = u;
			}
		}

		while(Q_n > 0&&keys[Q[Q_n-1]] < K) {
			ui u = Q[-- Q_n];

			ui pos = cstart[c_n];
			ids[pos ++] = u; vis[u] = 4;
			while(sv_next[u] != u) {
				u = sv_next[u];
				ids[pos ++] = u; vis[u] = 4;
			}
			cstart[++ c_n] = pos;
		}

		for(ui i = 0;i < Q_n;i ++) vis[Q[i]] = 0;
		if(Q_n == 0) break;
		//printf("finished one iteration\n\n");
	}
	for(ui i = old_cn;i < c_n;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) vis[ids[j]] = 0;

	return c_n;
}

void Graph::merge(Edge **graph_head, Edge *edges, ui u, ui v, ui *keys, ui *sv_next, ui *sv_last) {
	for(Edge *e = graph_head[v];e != nullptr;) {
		Edge *tmp = e->next;

		if(e->vertex == u) {
			-- keys[u];
			delete_edge(u, e->reverse, graph_head);
		}
		else {
			assert(e->reverse->vertex == v);
			e->reverse->vertex = u;

			add_edge(u, e, graph_head);
		}

		e = tmp;
	}
	graph_head[v] = nullptr;

	sv_next[sv_last[u]] = v;
	sv_last[u] = sv_last[v];
}

void Graph::delete_edge(ui u, Edge *e, Edge **graph_head) {
	if(e->pre == nullptr) {
		assert(graph_head[u] == e);
		e = e->next;
		if(e != nullptr) e->pre = nullptr;
		graph_head[u] = e;
	}
	else {
		assert(graph_head[u] != e);
		e->pre->next = e->next;
		if(e->next != nullptr) e->next->pre = e->pre;
	}
}

void Graph::add_edge(ui u, Edge *e, Edge **graph_head) {
	if(graph_head[u] != nullptr) graph_head[u]->pre = e;
	e->next = graph_head[u];
	graph_head[u] = e;
	e->pre = nullptr;
}

void Graph::remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *component_id, ui *pend, ui *Q, char *computed, ui *degree) {
	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) component_id[ids[j]] = i;

	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) {
		ui u = ids[j];
		for(ui k = pstart[u];k < pend[u];) {
			ui v = edges[k];
			assert(!computed[v]);
			if(component_id[v] != component_id[u]) swap(edges[k], edges[-- pend[u]]);
			else ++ k;
		}
	}

	ui Q_n = 0;
	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) {
		ui u = ids[j];
		degree[u] = pend[u] - pstart[u];
		if(degree[u] < K) {
			Q[Q_n ++] = u;
			computed[u] = 1;
		}
	}

	k_core_prune(K, Q, Q_n, computed, degree, pend) ;
}

void Graph::remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *pend, ui *Q, char *computed, ui *degree, UnionFind *UF) {
	ui Q_n = 0;
	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) {
		ui u = ids[j];
		ui ru = UF->UF_find(u);
		for(ui k = pstart[u];k < pend[u];) {
			if(UF->UF_find(edges[k]) != ru) swap(edges[k], edges[-- pend[u]]);
			else ++ k;
		}
		degree[u] = pend[u] - pstart[u];
		if(degree[u] < K) {
			Q[Q_n ++] = u;
			computed[u] = 1;
		}
	}
	k_core_prune(K, Q, Q_n, computed, degree, pend);
}

void Graph::print_kecc(ui K, ui c_n, ui *cstart, ui *ids) {
	ostringstream os;
	os<<dir<<"/"<<K<<"-ECC.txt";
	FILE *fout = Utility::open_file(os.str().c_str(), "w");

	for(ui i = 0;i < c_n;i ++) {
		sort(ids+cstart[i], ids+cstart[i+1]);
		for(ui j = cstart[i];j < cstart[i+1];j ++) fprintf(fout, "%u ", ids[j]);
		fprintf(fout, "\n");
	}

	fclose(fout);
}
