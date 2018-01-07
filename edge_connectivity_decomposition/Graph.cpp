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

void Graph::edge_connectivity_decomposition(bool print) {
#ifdef NDEBUG
	printf("*** ec_decomposition (Release): %s ***\n", dir.c_str());
#else
	printf("*** ec_decomposition (Debug): %s ***\n", dir.c_str());
#endif

	read_graph_binary();

	Timer timer;

	ui *Q = new ui[n];
	ui *active_component = new ui[n];
	ui max_core = core_decomposition(Q, active_component);

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	stack<pair<pair<ui,ui>, pair<ui,ui> > > working_stack; // (vertex, active_id of parent computation, L, H)

	for(ui i = 0;i < n;i ++) if(!vis[i]) {
		ui Q_n = 1;
		Q[0] = i; vis[i] = 1;
		for(ui j = 0;j < Q_n;j ++) {
			ui u = Q[j];
			for(ui k = pstart[u];k < pstart[u+1];k ++) if(!vis[edges[k]]) {
				vis[edges[k]] = 1;
				Q[Q_n ++] = edges[k];
			}
		}

		for(ui j = 0;j < Q_n;j ++) active_component[Q[j]] = i;
		working_stack.push(make_pair(make_pair(i,i), make_pair(2, max_core)));
	}

	UnionFind *UF = new UnionFind(n);
	UF->init(n);

	ui *representative = new ui[n]; // representative vertex of super vertex
	for(ui i = 0;i < n;i ++) representative[i] = i;

	ui *sv_next = new ui[n]; // next of super vertex
	ui *sv_last = new ui[n]; // last of super vertex
	for(ui i = 0;i < n;i ++) sv_next[i] = sv_last[i] = i;

	while(!working_stack.empty()) {
		ui u = working_stack.top().first.first, active_parent = working_stack.top().first.second;
		ui L = working_stack.top().second.first, H = working_stack.top().second.second;

		ui M = (L+H+1)/2;
		assert(L <= M <= H);

		if(L == M) working_stack.pop();
		else working_stack.top().second.second = M-1;


	}

	delete[] sv_next;
	delete[] sv_last;

	delete UF;
	delete[] representative;

	delete[] vis;
	delete[] Q;
	delete[] active_component;

	printf("\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
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

	printf("\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
}

void Graph::k_edge_connected_component_space(ui K, bool print) {
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

	char *computed = new char[n];
	memset(computed, 0, sizeof(char)*n);

	ui *Q = new ui[n], Q_n = 0;
	ui *active_component = new ui[n];
	// set active component id for vertices
	for(ui i = 0;i < n;i ++) if(!computed[i]) {
		Q[0] = i; Q_n = 1; computed[i] = 1;
		for(ui j = 0;j < Q_n;j ++) {
			ui u = Q[j]; active_component[u] = i;
			for(ui k = pstart[u];k < pstart[u+1];k ++) if(!computed[edges[k]]) {
				Q[Q_n ++] = edges[k];
				computed[edges[k]] = 1;
			}
		}
	}

	UnionFind *UF = new UnionFind(n);
	UF->init(n);

	ui *representative = new ui[n];
	for(ui i = 0;i < n;i ++) representative[i] = i;

	ui *sv_next = new ui[n];
	ui *sv_last = new ui[n];
	ui *adj_next = new ui[n];
	for(ui i = 0;i < n;i ++) sv_next[i] = sv_last[i] = adj_next[i] = i;

	ui *degree = new ui[n];
	// k-core-based pruning
	memset(computed, 0, sizeof(char)*n);
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1] - pstart[i];
		if(degree[i] < K) {
			Q[Q_n ++] = i;
			computed[i] = 1;
		}
	}
	//k_core_prune(Q, Q_n, computed, degree, active_component);

	ui *ids = new ui[n];
	ui *cstart = new ui[n+1];
	ui c_n = 0;
	cstart[0] = 0;

	ListLinearHeap *heap = new ListLinearHeap(n, K);

	for(ui i = 0;i < n;) {
		if(computed[i]) {
			++ i;
			continue;
		}

		ui new_cn = c_n;
		while(true) {
			if(adj_next[i] == i&&active_component[representative[UF->UF_find(edges[pstart[i]])]] != active_component[i]) {
				ui pos = cstart[new_cn];
				ui u = i; ids[pos ++] = u;
				while(sv_next[u] != u) {
					u = sv_next[u];
					ids[pos ++] = u;
				}
				cstart[++ new_cn] = pos;
				break;
			}

			heap->init(0, K, nullptr, nullptr);
			heap->insert(i, 0);

			ui u, key;
			while(heap->pop_max(u, key)) {

			}
		}
		if(new_cn == c_n + 1) {
			for(ui j = cstart[new_cn-1];j < cstart[new_cn];j ++) computed[ids[j]] = 1;
			if(cstart[new_cn] > cstart[new_cn-1] + 1) ++ c_n;
		}
	}

	delete[] ids;
	delete[] cstart;

	delete[] representative;
	delete UF;
	delete[] sv_next;
	delete[] sv_last;
	delete[] adj_next;

	delete[] active_component;
	delete[] computed;
	delete[] degree;
	delete[] Q;

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

ui Graph::decomposition(ui s, ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, Edge **graph_head, Edge *graph_edges, ListLinearHeap *heap, ui *sv_next, ui *sv_last) {
	//decomposition
	while(graph_head[s] != nullptr) {
		heap->init(0, K, nullptr, nullptr);
		heap->insert(s, 0);
		ui u, key, Q_n = 0;
		//printf("start an iteration\n");
		while(heap->pop_max(u, key)) {
			Q[Q_n ++] = u; vis[u] = 1; //vis[u] = 1 means u is in Q, and vis[u] = 2 means u is in heap
										//vis[u] == 3 means u is in Q between Q_n and new_Qn
			keys[u] = key;
			//printf(" [%u,%u]\n", u, key);

			ui new_Qn = Q_n;
			for(ui j = Q_n-1;j < new_Qn;j ++) {
				ui v = Q[j];
				//printf(" %d\n", v);
				//if(v == 0) printf("graph_head[0]: %u\n", graph_head[0]);

				for(Edge *e = graph_head[v];e != nullptr;e = e->next) {
					//if(v == 0) printf("graph_head[0]: %u\n", graph_head[0]);

					//printf("before\n");
					ui w = e->vertex;
					//printf("after\n");
					//if(v == 0) printf("%u is an neighbor of %u\n", w, v);

					if(vis[e->vertex] != 1) {
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
				}

				if(v == u) continue;

				// contract u and v
				vis[v] = 0;
				keys[u] += keys[v];
				//printf("merged %u and %u\n", u, v);
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
		//printf("finished one iteration\n\n");
	}
	return c_n;
	//printf("finished decomposition\n");
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
