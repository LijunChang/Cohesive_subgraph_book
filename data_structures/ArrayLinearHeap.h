/*
 * ArrayLinearHeap.h
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef DATA_STRUCTURES_ARRAYLINEARHEAP_H_
#define DATA_STRUCTURES_ARRAYLINEARHEAP_H_

#include "../utilities/Defines.h"

class ArrayLinearHeap {
private:
	ui n; // number vertices
	ui key_cap; // the maximum allowed key value

	ui max_key; // possible max key
	ui min_key; // possible min key

	ui *keys; // key of vertices

	ui *heads; // head of array of all vertices with a specific weight
	ui *ids; // ids of vertices in non-decreasing key order
	ui *rids; // reverse of id, i.e., the position of vertex in the ids array

public:
	ArrayLinearHeap(ui _n, ui _key_cap) {
		n = _n;
		key_cap = _key_cap;

		max_key = 0; min_key = key_cap;

		heads = ids = rids = keys = nullptr;
	}

	~ArrayLinearHeap() {
		if(heads != nullptr) {
			delete[] heads;
			heads = nullptr;
		}
		if(ids != nullptr) {
			delete[] ids;
			ids = nullptr;
		}
		if(rids != nullptr) {
			delete[] rids;
			rids = nullptr;
		}
		if(keys != nullptr) {
			delete[] keys;
			keys = nullptr;
		}
	}

	void init(ui _n, ui _key_cap, ui *_ids, ui *_keys) {
		if(keys == nullptr) keys = new ui[n];
		if(ids == nullptr) ids = new ui[n];
		if(rids == nullptr) rids = new ui[n];
		if(heads == nullptr) heads = new ui[key_cap+2];

		max_key = 0; min_key = _key_cap;

		ui *cnt = heads;
		memset(cnt, 0, sizeof(ui)*(_key_cap+1));
		for(ui i = 0;i < _n;i ++) {
			ui key = _keys[i];
			keys[_ids[i]] = key;
			assert(key <= _key_cap);
			++ cnt[key];

			if(key > max_key) max_key = key;
			if(key < min_key) min_key = key;
		}
		for(ui i = 1;i <= max_key;i ++) cnt[i] += cnt[i-1];
		for(ui i = 0;i < _n;i ++) rids[_ids[i]] = -- cnt[keys[i]];
		for(ui i = 0;i < _n;i ++) ids[rids[_ids[i]]] = _ids[i];

		for(ui i = 0, j = 0;i <= max_key + 1;i ++) {
			while(j < _n&&keys[ids[j]] < i) ++ j;
			heads[i] = j;
		}
	}

	ui get_n() { return n; }
	ui get_key_cap() { return key_cap; }
	ui get_key(ui id) { return keys[id]; }

	bool empty() {
		while(min_key <= max_key&&heads[min_key] >= heads[min_key+1]) ++ min_key;
		while(min_key <= max_key&&heads[max_key] >= heads[max_key+1]) -- max_key;

		return min_key > max_key;
	}

	// get the (id,key) pair with the maximum key value; return true if success, return false otherwise
	bool get_max(ui &id, ui &key) {
		if(empty()) return false;

		id = ids[heads[max_key+1] - 1];
		key = max_key;
		assert(keys[id] == key);

		return true;
	}

	// pop the (id,key) pair with the maximum key value; return true if success, return false otherwise
	bool pop_max(ui &id, ui &key) {
		if(empty()) return false;

		id = ids[-- heads[max_key+1]];
		key = max_key;
		assert(keys[id] == key);

		return true;
	}

	// get the (id,key) pair with the minimum key value; return true if success, return false otherwise
	bool get_min(ui &id, ui &key) {
		if(empty()) return false;

		id = ids[heads[min_key]];
		key = min_key;
		assert(keys[id] == key);

		return true;
	}

	// pop the (id,key) pair with the minimum key value; return true if success, return false otherwise
	bool pop_min(ui &id, ui &key) {
		if(empty()) return false;

		id = ids[heads[min_key] ++];
		key = min_key;
		assert(keys[id] == key);

		return true;
	}

	// increment the key of vertex id by inc
	ui increment(ui id, ui inc) {
		assert(keys[id] + inc <= key_cap);

		while(inc --) {
			ui &key = keys[id];
			ui pos1 = heads[key+1]-1, pos2 = rids[id];
			std::swap(ids[pos1], ids[pos2]);
			rids[ids[pos1]] = pos1; rids[ids[pos2]] = pos2;

			if(max_key == key) {
				++ max_key;
				heads[max_key+1] = heads[max_key];
			}

			++ key;
			-- heads[key];
		}
		return keys[id];
	}

	// decrement the key of vertex id by dec
	ui decrement(ui id, ui dec) {
		assert(keys[id] >= dec);

		while(dec --) {
			ui &key = keys[id];
			ui pos1 = heads[key], pos2 = rids[id];
			std::swap(ids[pos1], ids[pos2]);
			rids[ids[pos1]] = pos1; rids[ids[pos2]] = pos2;

			if(min_key == key) {
				-- min_key;
				heads[min_key] = heads[min_key+1];
			}

			++ heads[key];
			-- key;
		}
		return keys[id];
	}
};

#endif /* DATA_STRUCTURES_ARRAYLINEARHEAP_H_ */
