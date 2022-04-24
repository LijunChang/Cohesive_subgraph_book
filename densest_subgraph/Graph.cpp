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

	original_n = 0;
	component_id = nullptr;

	pstart = nullptr;
	pend = nullptr;
	edges = nullptr;

	res.clear();

	revers = nullptr;
	capacity_pre = nullptr;
	capacity = nullptr;
	level = nullptr;
	now = nullptr;
	que = nullptr;

	total_flow_pre = lambda_pre = 0;
}

Graph::~Graph() {
	if(component_id != nullptr) {
		delete[] component_id;
		component_id = nullptr;
	}
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(revers != nullptr) {
		delete[] revers;
		revers = nullptr;
	}
	if(capacity_pre != nullptr) {
		delete[] capacity_pre;
		capacity_pre = nullptr;
	}
	if(capacity != nullptr) {
		delete[] capacity;
		capacity = nullptr;
	}
	if(level != nullptr) {
		delete[] level;
		level = nullptr;
	}
	if(now != nullptr) {
		delete[] now;
		now = nullptr;
	}
	if(que != nullptr) {
		delete[] que;
		que = nullptr;
	}
}

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
	pstart = new ui[n+3];
	if(pend != nullptr) delete[] pend;
	pend = new ui[n+3];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m+4*n];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
		pstart[i+1] = pend[i] = pstart[i] + degree[i];
	}

	fclose(f);

	delete[] degree;

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

void Graph::read_graph(const char *input_file) {
	printf("# Start reading graph %s\n", input_file);
	FILE *f = Utility::open_file(input_file, "r");

	fscanf(f, "%u%u", &n, &m);
	++ n; m *= 2;
	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	vector<pair<ui,ui> > vp;
	for(ui i = 0;i < m/2;i ++) {
		ui a, b;
		fscanf(f, "%u%u", &a, &b);
		vp.pb(mp(a,b));
		vp.pb(mp(b,a));
	}
	sort(vp.begin(), vp.end());

	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[n+3];
	if(pend != nullptr) delete[] pend;
	pend = new ui[n+3];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m+4*n];

	pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		pend[i] = pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) edges[pend[i] ++] = vp[idx ++].second;
		pstart[i+1] = pend[i];
	}

	fclose(f);

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

void Graph::densest_subgraph() {
	original_n = n;

	if(m == 0) {
		res.pb(1);
		return ;
	}

	Timer timer;
	ui *peel_sequence = new ui[n+2];
	ui *core = new ui[n+2];
	ui *degree = new ui[n+2];
	ui *rid = new ui[n+2];
	vector<ui> ids;

	ui best_n, best_m;
	ui max_core = core_decomposition(peel_sequence, core, n, pstart, pend, edges);
	greedy_and_reduce(best_n, best_m, peel_sequence, core, ids, degree, rid);

	// check 2m*n*(n-1) fits in lint
	lint tmp = 2*(lint)m*(lint)n*(lint)(n-1);
	if(tmp/2/m/n != n-1) {
		printf("2m*n*(n-1) does not fit in lint!\n");
		delete[] degree;
		delete[] peel_sequence;
		delete[] core;
		delete[] rid;

		return ;
	}

	assert(revers == nullptr);
	revers = new ui[2*m+4*n];
	construct_lambda_graph();

	assert(capacity_pre == nullptr&&capacity == nullptr);
	capacity_pre = new lint[2*m+4*n];
	capacity = new lint[2*m+4*n];
	for(ui i = 0;i < n;i ++) {
		for(ui j = pstart[i];j < pend[i]-2;j ++) capacity_pre[j] = (lint)n*(lint)(n-1);
		capacity_pre[pend[i]-2] = 0; capacity_pre[pstart[n]+i] = (lint)degree[i]*(lint)n*(lint)(n-1);
		capacity_pre[pend[i]-1] = 0; capacity_pre[pstart[n+1]+i] = 0;
	}
	total_flow_pre = 0;
	lambda_pre = 0;

	assert(level == nullptr); level = new ui[n+2];
	assert(que == nullptr); que = new ui[n+2];
	assert(now == nullptr); now = new ui[n+2];
	char *vis = new char[n+2];

	ui x = best_m/best_n;
	lint low = (lint)(x)*(lint)(n)*(n-1) + ((lint)(best_m-x*best_n)*(lint)(n)*(lint)(n-1)+best_n-1)/(lint)(best_n);
	lint high = (lint)max_core*(lint)(n)*(n-1);

	lint total_flow= min_cut(1, low);
	if(total_flow < 2*m*(lint)n*(lint)(n-1)) {
		save_S(res, vis);

		for(ui i = 0;i < pend[n+1];i ++) capacity_pre[i] = capacity[i];
		total_flow_pre = total_flow; lambda_pre = low;

		++ low;
	}
	else low = high+1;

	ui searches = 1;

	while (low < high) {
		lint mid = low + (high - low) / 2;
		//printf("checked %lld\n", mid);
		total_flow = min_cut(1, mid);
		++ searches;
		if (total_flow < 2 * m * (lint) n * (lint) (n - 1)) {
			low = mid + 1;
			save_S(res, vis);

			for (ui i = 0; i < pend[n + 1]; i++) capacity_pre[i] = capacity[i];
			total_flow_pre = total_flow; lambda_pre = mid;
		}
		else high = mid - 1;
	}
	sort(res.begin(), res.end());

	memset(vis, 0, sizeof(char)*(n+2));
	for(ui i = 0;i < res.size();i ++) vis[res[i]] = 1;
	best_n = res.size(); best_m = 0;
	for(ui i = 0;i < res.size();i ++) for(ui j = pstart[res[i]];j < pend[res[i]];j ++) if(vis[edges[j]]) ++ best_m;
	best_m /= 2;

	printf("Densest subgraph: n = %u, m = %s, number of searches: %u\n", best_n, Utility::integer_to_string(best_m).c_str(), searches);

	if(res.size() < n) {
#ifndef NDEBUG
		for(ui i = 1;i < res.size();i ++) assert(res[i] > res[i-1]);
		assert(n == ids.size());
		for(ui i = 1;i < n;i ++) assert(ids[i] > ids[i-1]);
#endif

		for(ui i = 0;i < n+2;i ++) rid[i] = n;
		for(ui i = 0;i < res.size();i ++) rid[res[i]] = i;

		n = res.size();
		m = 0;
		for(ui i = 0;i < n;i ++) {
			ui old_start = pstart[res[i]], old_end = pend[res[i]];
			if(i) pstart[i] = pstart[i-1] + degree[i-1] + 2;
			else pstart[i] = 0;
			pend[i] = pstart[i];

			degree[i] = 0;
			for(ui j = old_start;j < old_end;j ++) if(rid[edges[j]] < n) {
				++ degree[i];
				if(rid[edges[j]] > i) edges[pend[i] ++] = rid[edges[j]];
			}
			m += degree[i];
		}
		pstart[n] = pstart[n-1] + degree[n-1] + 2;
		m /= 2;

		printf("reduced for construct index\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m).c_str());

		construct_lambda_graph();

		for(ui i = 0;i < n;i ++) res[i] = ids[res[i]];
		ids = res;
		for(ui i = 0;i < n;i ++) res[i] = i;
	}

	for(ui i = 0;i < n;i ++) {
		assert(edges[pend[i]-2] == n&&edges[pend[i]-1] == n+1);
		for(ui j = pstart[i];j < pend[i]-2;j ++) capacity[j] = (lint)n;
		capacity[pend[i]-2] = 0; capacity[pstart[n]+i] = (lint)degree[i]*(lint)n;
		capacity[pend[i]-1] = 2*m; capacity[pstart[n+1]+i] = 0;
	}

	total_flow = min_cut(0, 1);
	assert(total_flow == 2*m*(lint)n);

	Timer t2;

	for(ui i = 0;i < n;i ++) {
		ui old_end = pend[i]-2, end = pstart[i];
		assert(edges[old_end] == n&&edges[old_end+1] == n+1);
		for(ui j = pstart[i];j < old_end;j ++) if(capacity[j] > 0) edges[end ++] = edges[j];
		pend[i] = end;
	}

	ui scc_n;
	ui *scc_starts = new ui[res.size()+1];
	ui *scc_ids = new ui[res.size()];
	compute_SCCs(scc_n, scc_starts, scc_ids, n, pstart, pend, edges);

	ui *component_id = degree;
	for(ui i = 0;i < scc_n;i ++) for(ui j = scc_starts[i];j < scc_starts[i+1];j ++) component_id[scc_ids[j]] = i;

	ui minimal_densest_n = 0;
	for(ui i = 0;i < scc_n;i ++) {
		ui ok = 1;
		for(ui j = scc_starts[i];j < scc_starts[i+1]&&ok;j ++) {
			ui u = scc_ids[j];
			for(ui k = pstart[u];k < pend[u];k ++) {
				assert(edges[k] >= 0&&edges[k] < n);
				if(component_id[edges[k]] != i) {
					ok = 0;
					break;
				}
			}
		}
		if(ok) ++ minimal_densest_n;
	}

	printf("Total %u minimal densest subgraphs\n", minimal_densest_n);
	if(scc_n == 1) printf("Minimal is the same as the maximal\n");

	printf("\tIndex construction time: %s\n", Utility::integer_to_string(t2.elapsed()).c_str());

	for(ui i = 0;i < res.size();i ++) res[i] = ids[res[i]];

#ifndef NDEBUG
	for(ui i = 1;i < res.size();i ++) assert(res[i] > res[i-1]);
#endif

	delete[] scc_starts;
	delete[] scc_ids;
	delete[] degree;
	delete[] core;
	delete[] peel_sequence;
	delete[] rid;
	delete[] vis;

	printf("\tTotal processing time excluding I/O: %s\n\n", Utility::integer_to_string(timer.elapsed()).c_str());
}

void Graph::print_densest_subgraph() {
	FILE *f = Utility::open_file((dir + string("/densest_subgraph.txt")).c_str(), "w");
	fprintf(f, "The maximal densest subgraph contains %lu vertices\n", res.size());
	for(ui i = 0;i < res.size();i ++) fprintf(f, "%u\n", res[i]);
	fclose(f);
}

void Graph::print_densest_subgraph(const char *output_file) {
	FILE *f = Utility::open_file(output_file, "w");
	fprintf(f, "%lu\n", res.size());
	for(ui i = 0;i < res.size();i ++) fprintf(f, "%u\n", res[i]);
	fclose(f);
}

// private member functions

// core decomposition:
//		"peel_sequence" stores the vertices in degenerarcy order
//		"core[u]" stores the core value of vertex u
//		"n, pstart, pend, edges" are used to represent the graph
// return the maximum core value
ui Graph::core_decomposition(ui *peel_sequence, ui *core, const ui n, const ui *pstart, const ui *pend, const ui *edges) {
#ifndef NDEBUG
	printf("Start core decomposition\n");
#endif
	Timer timer;
	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pend[i]-pstart[i];

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
		for(ui j = pstart[u];j < pend[u];j ++) if(rid[edges[j]] > i) {
			ui v = edges[j];
			ui pos1 = degree_start[degree[v]], pos2 = rid[v];
			std::swap(id[pos1], id[pos2]);
			rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
			++ degree_start[degree[v]];
			-- degree[v];
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) {
		ui cnt = 0;
		ui u = peel_sequence[i];
		for(ui j = pstart[u];j < pend[u];j ++) if(rid[edges[j]] > i) ++ cnt;
		assert(cnt == degree[u]);
	}
#endif

	delete[] degree;
	delete[] degree_start;
	delete[] rid;

	printf("core decomposition time: %s (microseconds)\n", Utility::integer_to_string(timer.elapsed()).c_str());

#ifndef NDEBUG
	printf("Finished core decomposition\n");
#endif
	return max_core;
}

// compute SCCs of a directed graph, and represent the SCCs by scc_starts and scc_ids
// Tarjan's one-pass DFS algorithm
// the SCC IDs are assigned in reverse topological sort order
void Graph::compute_SCCs(ui &scc_n, ui *scc_starts, ui *scc_ids, const ui n, const ui *pstart, const ui *pend, const ui *edges) {
	ui *index = new ui[n];
	ui *low = new ui[n];
	ui *pos = new ui[n];
	ui *stk = new ui[n];
	ui *parent = new ui[n];
	char *in_stk = new char[n];

	for(ui i = 0;i < n;i ++) index[i] = n;
	for(ui i = 0;i < n;i ++) pos[i] = pstart[i];
	memset(in_stk, 0, sizeof(char)*n);

	ui cnt = 0; scc_n = 0; scc_starts[0] = 0;
	for(ui i = 0;i < n;i ++) if(index[i] == n) {
		index[i] = cnt;
		low[i] = cnt;
		++ cnt;
		ui stk_n = 1;
		stk[0] = i;
		in_stk[i] = 1;

		ui u = i;
		while(stk_n) {
			while(pos[u] < pend[u]) {
				ui v = edges[pos[u]];
				if(index[v] == n) {
					index[v] = cnt;
					low[v] = cnt;
					++ cnt;
					stk[stk_n ++] = v;
					in_stk[v] = 1;
					parent[v] = u;
					u = v;
					break;
				}
				if(in_stk[v]) low[u] = min(low[u], low[v]);
				++ pos[u];
			}

			if(pos[u] >= pend[u]) {
				if(index[u] == low[u]) {
					++ scc_n; scc_starts[scc_n] = scc_starts[scc_n-1];
					while(stk_n) {
						ui v = stk[-- stk_n];
						scc_ids[scc_starts[scc_n] ++] = v;
						in_stk[v] = 0;
						if(v == u) break;
					}
				}

				u = parent[u];
			}
		}
	}

	delete[] in_stk;
	delete[] parent;
	delete[] index;
	delete[] low;
	delete[] pos;
	delete[] stk;
}

// the resulting graph only store the neighbors with higher IDs than itself
void Graph::greedy_and_reduce(ui &best_n, ui &best_m, ui *peel_sequence, ui *core, vector<ui> &ids, ui *degree, ui *rid) {
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;

	// greedy 2-approximate densest subgraph
	m /= 2; best_n = n; best_m = m;
	ui best_idx = 0, current_m = m;
	for(ui i = 0;i < n;i ++) {
		if((lint)current_m*(lint)best_n > (lint)best_m*(lint)(n-i)) {
			best_n = n-i;
			best_m = current_m;
			best_idx = i;
		}
		for(ui j = pstart[peel_sequence[i]];j < pend[peel_sequence[i]];j ++) if(rid[edges[j]] > i) -- current_m;
	}
	res.clear();
	for(ui i = best_idx;i < n;i ++) res.pb(peel_sequence[i]);

	// reduce to [best_m/best_n]-core
	ui start_idx = 0, threshold = (best_m+best_n-1)/best_n;
	while(start_idx < n&&core[peel_sequence[start_idx]] < threshold) ++ start_idx;
	//printf("start_idx: %u, best_idx: %u, n: %u\n", start_idx, best_idx, n);
	assert(start_idx <= best_idx);
	for(ui i = start_idx;i < n;i ++) peel_sequence[i-start_idx] = peel_sequence[i];
	for(ui i = 0;i < n;i ++) rid[i] = n;
	n -= start_idx;

	// relabel vertex IDs
	ids.clear(); ids.reserve(n);
	for(ui i = 0;i < n;i ++) ids.pb(peel_sequence[i]);
	sort(ids.begin(), ids.end());
	for(ui i = 0;i < n;i ++) rid[ids[i]] = i;
	for(ui i = 0;i < n;i ++) {
		peel_sequence[i] = rid[peel_sequence[i]];
		core[i] = core[ids[i]];
	}
	for(ui i = 0;i < res.size();i ++) res[i] = rid[res[i]];

	// get the degree of vertices in the reduced graph
	m = 0;
	for(ui i = 0;i < n;i ++) {
		degree[i] = 0;
		for(ui j = pstart[ids[i]];j < pend[ids[i]];j ++) if(rid[edges[j]] < n) ++ degree[i];
		m += degree[i];
	}
	m /= 2;

	printf("reduced (%u-core): n = %s; m = %s (undirected)\n", threshold, Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m).c_str());

	// get higher ID neighbors
	ui *tmp = new ui[2*m+4*n];
	for(ui i = 0;i < n;i ++) {
		ui old_start = pstart[ids[i]], old_end = pend[ids[i]];
		if(i) pstart[i] = pstart[i-1] + degree[i-1] + 2;
		else pstart[i] = 0;
		pend[i] = pstart[i];

		for(ui j = old_start;j < old_end;j ++) if(rid[edges[j]] < n&&rid[edges[j]] > i) {
			tmp[pend[i] ++] = rid[edges[j]];
		}
	}
	pstart[n] = pstart[n-1] + degree[n-1] + 2;
	delete[] edges;
	edges = tmp;
	tmp = nullptr;
}

void Graph::construct_lambda_graph() {
	for(ui i = 0;i < n;i ++) {
		for(ui j = pstart[i];j < pend[i]&&edges[j] > i;j ++) {
			ui v = edges[j];
			revers[j] = pend[v]; revers[pend[v]] = j;
			edges[pend[v] ++] = i;
		}
	}
	pend[n] = pstart[n];
	for(ui i = 0;i < n;i ++) {
		revers[pend[i]] = pend[n]; revers[pend[n]] = pend[i];
		edges[pend[i] ++] = n; edges[pend[n] ++] = i;
	}
	pend[n+1] = pstart[n+1] = pend[n];
	for(ui i = 0;i < n;i ++) {
		revers[pend[i]] = pend[n+1]; revers[pend[n+1]] = pend[i];
		edges[pend[i] ++] = n+1; edges[pend[n+1] ++] = i;
	}
	pstart[n+2] = pend[n+1];
}

bool Graph::bfs() {
	memset(level, 0, sizeof(ui)*(n+2));
	for(ui i = 0;i < n+2;i ++) now[i] = pstart[i];

	ui queue_n = 1;
	que[0] = n;
	level[n] = 1;

	for(ui i = 0;i < queue_n;i ++) {
		ui u = que[i];
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];
			if(capacity[j]&&level[v] == 0) {
				level[v] = level[u] + 1;
				que[queue_n ++] = v;
				if(v == n+1) return true;
			}
		}
	}

	return false;
}

lint Graph::dinic(ui u, lint flow) {
	if(u == n+1) return flow;
	lint pre_flow = flow;
	for(;now[u] < pend[u]&&flow;now[u] ++) {
		ui idx = now[u];
		ui v = edges[idx];
		if(capacity[idx] == 0||level[v] != level[u] + 1 ) continue;
		lint tmp = dinic(v, min(capacity[idx], flow));
		capacity[idx] -= tmp;
		capacity[revers[idx]] += tmp;
		flow -= tmp;
	}
	if(pre_flow == flow) level[u] = n+10;
	return pre_flow - flow;
}

lint Graph::min_cut(char initialize, lint lambda) {
	lint total_flow = 0;
	if(initialize) {
		for(ui i = 0;i < pend[n+1];i ++) capacity[i] = capacity_pre[i];
		for(ui i = 0;i < n;i ++) capacity[pend[i]-1] += 2*(lambda - lambda_pre);
		total_flow += total_flow_pre;
	}

	while(bfs()) total_flow += dinic(n, LINT_INF);

	return total_flow;
}

void Graph::save_S(vector<ui> &res, char *vis) {
	res.clear();
	memset(vis, 0, sizeof(char)*(n+2));
	for(ui i = pstart[n];i < pend[n];i ++) if(capacity[i]) {
		res.pb(edges[i]);
		vis[edges[i]] = 1;
	}
	for(ui i = 0;i < res.size();i ++) {
		ui u = res[i];
		for(ui j = pstart[u];j < pend[u]-2;j ++) if(capacity[j]&&!vis[edges[j]]) {
			res.pb(edges[j]);
			vis[edges[j]] = 1;
		}
	}
	//printf("save size: %lu\n", res.size());
}
