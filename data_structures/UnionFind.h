/*
 * UnionFind.h
 *
 *  Created on: 6Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef DATA_STRUCTURES_UNIONFIND_H_
#define DATA_STRUCTURES_UNIONFIND_H_

#include "../utilities/Defines.h"

class UnionFind {
private:
	ui n; // number of elements
	ui *parent;
	ui *rank;

public:
	UnionFind(ui _n) {
		n = _n;
		parent = rank = nullptr;
	}
	~UnionFind() {
		if(parent != nullptr) {
			delete[] parent;
			parent = nullptr;
		}
		if(rank != nullptr) {
			delete[] rank;
			rank = nullptr;
		}
	}

	void init(ui _n = 0) {
		if(_n == 0) _n = n;
		assert(_n <= n);

		if(parent == nullptr) parent = new ui[n];
		if(rank == nullptr) rank = new ui[n];

		for(ui i = 0;i < _n;i ++) {
			parent[i] = i;
			rank[i] = 0;
		}
	}

	void init(ui *ids, ui _n) {
		assert(_n <= n);

		if(parent == nullptr) parent = new ui[n];
		if(rank == nullptr) rank = new ui[n];

		for(ui i = 0;i < _n;i ++) {
			ui u = ids[i];
			parent[u] = u;
			rank[u] = 0;
		}
	}

	ui UF_find(ui u) {
		ui res = u;
		while(parent[res] != res) res = parent[res];
		while(parent[u] != res) {
			ui tmp = parent[u];
			parent[u] = res;
			u = tmp;
		}
		return res;
	}

	// return the new root of the merged tree
	ui UF_union(ui u, ui v) {
		ui tu = UF_find(u);
		ui tv = UF_find(v);

		if(tu == tv) return tu;

		ui res;
		if(rank[tu] > rank[tv]) {
			res = tu;
			parent[tv] = tu;
		}
		else {
			res = tv;
			parent[tu] = tv;
			if(rank[tu] == rank[tv]) ++ rank[tv];
		}
		return res;
	}
};

#endif /* DATA_STRUCTURES_UNIONFIND_H_ */
