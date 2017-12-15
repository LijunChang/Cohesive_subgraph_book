/*
 * ListLinearHeap.h
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef DATA_STRUCTURES_LISTLINEARHEAP_H_
#define DATA_STRUCTURES_LISTLINEARHEAP_H_

#include "../utilities/Defines.h"

class ListLinearHeap {
private:
	ui n; // number vertices
	ui key_cap; // the maximum allowed key value

	ui max_key; // possible max key
	ui min_key; // possible min key

	ui *keys; // keys of vertices
			  // keys[i] > key_cap if vertex i is not in the data structure

	ui *heads; // head of doubly-linked list for a specific weight
	ui *pres; // pre for doubly-linked list
	ui *nexts; // next for doubly-linked list

public:
	ListLinearHeap(ui _n, ui _key_cap) {
		n = _n;
		key_cap = _key_cap;

		min_key = key_cap;
		max_key = 0;

		heads = keys = pres = nexts = nullptr;
	}
	~ListLinearHeap() {
		if(heads != nullptr) {
			delete[] heads;
			heads = nullptr;
		}
		if(pres != nullptr) {
			delete[] pres;
			pres = nullptr;
		}
		if(nexts != nullptr) {
			delete[] nexts;
			nexts = nullptr;
		}
		if(keys != nullptr) {
			delete[] keys;
			keys = nullptr;
		}
	}

	// initialize the data structure by (id, key) pairs
	// _n is the number of pairs, _key_cap is the maximum possible key value
	void init(ui _n, ui _key_cap, ui *_ids, ui *_keys) {
		if(keys == nullptr) keys = new ui[n];
		if(pres == nullptr) pres = new ui[n];
		if(nexts == nullptr) nexts = new ui[n];
		if(heads == nullptr) heads = new ui[key_cap+1];
		assert(_key_cap <= key_cap);
		min_key = _key_cap; max_key = 0;
		for(ui i = 0;i <= _key_cap;i ++) heads[i] = n;

		for(ui i = 0;i < _n;i ++) insert(_ids[i], _keys[i]);
	}

	// insert (id, key) pair into the data structure
	void insert(ui id, ui key) {
		assert(id < n); assert(key <= key_cap);
		//assert(keys[id] > key_cap);

		keys[id] = key; pres[id] = n; nexts[id] = heads[key];
		if(heads[key] != n) pres[heads[key]] = id;
		heads[key] = id;

		if(key < min_key) min_key = key;
		if(key > max_key) max_key = key;
	}

	// remove a vertex from the data structure
	ui remove(ui id) {
		assert(keys[id] <= max_key);
		if(pres[id] == n) {
			assert(heads[keys[id]] == id);
			heads[keys[id]] = nexts[id];
			if(nexts[id] != n) pres[nexts[id]] = n;
		}
		else {
			ui pid = pres[id];
			nexts[pid] = nexts[id];
			if(nexts[id] != n) pres[nexts[id]] = pid;
		}

		ui ret = keys[id];
		keys[id] = key_cap+1;
		return ret;
	}

	ui get_n() { return n; }
	ui get_key_cap() { return key_cap; }

	// get the key value of a vertex; return true if the vertex is in the data structure, return false otherwise
	bool get_key(ui id, ui &key) {
		key = keys[id];
		return key <= key_cap;
	}

	void get_ids(std::vector<ui> &ids) {
		ids.clear();
		tighten();
		for(ui i = min_key;i <= max_key;i ++) {
			for(ui id = heads[i];id != n;id = nexts[id]) {
				ids.pb(id);
			}
		}
	}

	void get_ids_keys(std::vector<ui> &ids, std::vector<ui> &_keys) {
		ids.clear(); _keys.clear();
		tighten();
		for(ui i = min_key;i <= max_key;i ++) {
			for(ui id = heads[i];id != n;id = nexts[id]) {
				ids.pb(id); _keys.pb(id);
			}
		}
	}

	bool empty() {
		tighten();
		return min_key > max_key;
	}

	ui size() {
		tighten();
		ui res = 0;
		for(ui i = min_key;i <= max_key;i ++) for(ui id = heads[i];id != n;id = nexts[id]) ++ res;
		return res;
	}

	// get the (id,key) pair with the maximum key value; return true if success, return false otherwise
	bool get_max(ui &id, ui &key) {
		if(empty()) return false;

		id = heads[max_key];
		key = max_key;
		assert(keys[id] == key);
		return true;
	}

	// pop the (id,key) pair with the maximum key value; return true if success, return false otherwise
	bool pop_max(ui &id, ui &key) {
		if(empty()) return false;

		id = heads[max_key];
		key = max_key;
		assert(keys[id] == key);

		keys[id] = key_cap + 1;
		heads[max_key] = nexts[id];
		if(heads[max_key] != n) pres[heads[max_key]] = n;
		return true;
	}

	// get the (id,key) pair with the minimum key value; return true if success, return false otherwise
	bool get_min(ui &id, ui &key) {
		if(empty()) return false;

		id = heads[min_key];
		key = min_key;
		assert(keys[id] == key);

		return true;
	}

	// pop the (id,key) pair with the minimum key value; return true if success, return false otherwise
	bool pop_min(ui &id, ui &key) {
		if(empty()) return false;

		id = heads[min_key];
		key = min_key;

		assert(keys[id] == key);

		keys[id] = key_cap + 1;
		heads[min_key] = nexts[id];
		if(heads[min_key] != n) pres[heads[min_key]] = n;
		return true;
	}

	// increment the key of vertex id by inc
	ui increment(ui id, ui inc) {
		assert(keys[id]+inc <= key_cap);

		ui new_key = keys[id] + inc;

		remove(id);
		insert(id, new_key);

		return new_key;
	}

	// decrement the key of vertex id by dec
	ui decrement(ui id, ui dec) {
		assert(keys[id] >= dec);

		ui new_key = keys[id] - dec;

		remove(id);
		insert(id, new_key);

		return new_key;
	}

private:
	void tighten() {
		while(min_key <= max_key&&heads[min_key] == n) ++ min_key;
		while(min_key <= max_key&&heads[max_key] == n) -- max_key;
	}
};

#endif /* DATA_STRUCTURES_LISTLINEARHEAP_H_ */
