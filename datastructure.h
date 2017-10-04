#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include "pch.h"
#include "mesh.h"

class GroupLink {
public:
	GroupLink(int ida, int idb, float worth, float volume) {
		if (ida > idb) {
			this->ida = idb;
			this->idb = ida;
		}
		else {
			this->ida = ida;
			this->idb = idb;
		}
		this->worth = worth;
		this->volume = volume;
	}
	int ida, idb;
	float worth;
	float volume;
	bool operator < (const GroupLink& tar) const
	{
		if (worth == tar.worth)return volume > tar.volume;
		return worth < tar.worth;
	}
	bool operator == (const GroupLink& tar) const
	{
		if (worth == tar.worth)return volume == tar.volume;
		return worth == tar.worth;
	}
	bool operator > (const GroupLink& tar) const
	{
		if (worth == tar.worth)return volume < tar.volume;
		return worth > tar.worth;
	}
};

class TouchInfo{
public:
	TouchInfo() {}
	TouchInfo(int e, Vector3 va, Vector3 vb, Vector3 dir, int knifeid) {
		this->e = e;
		this->va = va;
		this->vb = vb;
		this->dir = dir;
		this->knifeid = knifeid;
	}
	int e, knifeid;
	Vector3 va, vb, dir;
};

class Hash {
public:
	std::vector<unsigned int> hash;
	int hashlen = 0;
	const int partlen = 32;
	Hash() {
		hashlen = 0;
		hash = std::vector<unsigned int>();
	}
	Hash(std::vector<unsigned int> hash, int hashlen) {
		this->hashlen = hashlen;
		this->hash = hash;
	}
	void addHash(bool positive) {
		int hashlendiv = hashlen / partlen;
		int hashlenmod = hashlen % partlen;
		if (hashlenmod == 0)hash.push_back(0x0);
		if (positive) {
			int mask = 0x1 << hashlenmod;
			hash[hashlendiv] = hash[hashlendiv] | mask;
		}
		hashlen++;
	}
	int difNum(Hash tar) {
		int cnt = 0;
		int lasthashlen = hashlen;
		for (int i = 0; i < hash.size(); i++, lasthashlen -= partlen)
		{
			unsigned int hashiter = hash[i] ^ tar.hash[i];
			int curhashlen = lasthashlen < partlen ? lasthashlen : partlen;
			for (int j = 0; j < curhashlen; j++, hashiter /= 2) {
				if (hashiter % 2 == 0x1)cnt++;
			}
		}
		return cnt;
	}
	bool getBit(int i) {
		unsigned int mask = 0x1 << (i%partlen);
		return (hash[i / partlen] & mask) != 0;
	}
	std::vector<bool> getBits() {
		std::vector<bool> bits = std::vector<bool>(hashlen);
		int flag = 0;
		//std::vector<bool> bits;
		int lasthashlen = hashlen;
		for (int i = 0; i < hash.size(); i++, lasthashlen -= partlen) {
			unsigned int hashiter = hash[i];
			int curhashlen = lasthashlen < partlen ? lasthashlen : partlen;
			for (int j = 0; j < curhashlen; j++, hashiter /= 2) {
				//if (hashiter % 2 == 0x1)bits.push_back(true);
				//else bits.push_back(false);
				if (hashiter % 2 == 0x1)bits[flag] = true;
				else bits[flag] = false;
				flag++;
			}
		}
		return bits;
	}
	Hash assign(std::vector<bool> &bits) {
		hashlen = bits.size();
		hash = std::vector<unsigned int>(bits.size() / partlen + 1);
		for (int i = 0; i < hash.size(); i++)hash[i] = 0x0;
		unsigned int mask = 0x1;
		for (int i = 0; i < bits.size(); i++) {
			if (bits[i])hash[i / partlen] = hash[i / partlen] | mask;
			mask = (mask << 1) | (mask >> (partlen - 1));
		}
		return *this;
	}
	bool operator < (const Hash& tar) const
	{
		for (int i = hash.size()-1; i >=0 ; i--) {
			if (hash[i] == tar.hash[i])continue;
			return hash[i] < tar.hash[i];
			
		}
		return false;
	}
	bool operator == (const Hash& tar) const
	{
		for (int i = hash.size() - 1; i >= 0; i--) {
			if (hash[i] != tar.hash[i])return false;
		}
		return true;
	}
	bool operator > (const Hash& tar) const
	{
		for (int i = hash.size() - 1; i >= 0; i--) {
			if (hash[i] == tar.hash[i])continue;
			return hash[i] > tar.hash[i];
		}
		return false;
	}
	void operator=(const Hash& hash) {
		this->hash = hash.hash;
		hashlen = hash.hashlen;
	}
	/*
	Hash() {
	hashlen = 0;
	hashes = std::vector<unsigned int>();
	}
	Hash(std::vector<unsigned int> hashes, int hashlen) {
	this->hashes = hashes;
	this->hashlen = hashlen;
	}
	int hashlen = 0;
	std::vector<unsigned int> hashes;
	void addHash(bool positive) {
	int hashpos = hashlen % 32;
	int hashgroup = hashlen / 32;
	if (hashpos == 0) {
	hashes.push_back(0);
	}
	if (positive) {
	int mask = 0x1 << hashpos;
	hashes[hashgroup] = hashes[hashgroup] | mask;
	}
	hashlen++;
	}
	*/
};

class Piece
{
public:
	Piece() {};
	Piece(Mesh mesh) {
		this->mesh = mesh;
	};
	Piece(Mesh mesh, Hash hash) {
		this->mesh = mesh;
		this->hash = hash;
	};
	Mesh mesh;
	Hash hash;
	int id;
	std::vector<TouchInfo> touchinfos;
	std::vector<Vector3> voxels;
	std::vector<int> voxelsi;
	bool exist = false;
	float volume = 0;
	bool isNei(Piece tar) {
		return hash.difNum(tar.hash) == 1;
	}
};

class Group
{
public:
	Group() {};
	//std::vector<Piece> pieces;
	std::vector<int> pieces;
	std::vector<Vector3> ava;
	int avanum=0;
	const float plus = 15.0f;
	int id;
	std::set<int> neighbor;//id
	float volume = 0;
	void initAva();
	Mesh getMesh(std::vector<Piece> & pieces);
	void operator=(const Group& group) {
		pieces = group.pieces;
		ava = group.ava;
		id = group.id;
		avanum = group.avanum;
		volume = group.volume;
	}
};

class Capsule {
public:
	Capsule(Vector3 p1, Vector3 p2, float r) {
		this->p1 = p1;
		this->p2 = p2;
		this->r = r;
	}
	Vector3 p1;
	Vector3 p2;
	float r;
	bool collide(Vector3 center, float size, bool &between) {
		Vector3 vec = (p2 - p1);
		float ll = vec.dot(vec);
		between = false;
		float dot = vec.normalize().dot(center - p1);//vec normalized
		if (dot < 0) {
			Vector3 d = (center - p1);
			return d.dot(d) < (r + size) * (r + size);
		}
		else if (dot*dot < ll) {
			between = true;
			Vector3 d = (center - p1);
			float h = vec.dot(d);
			return d.dot(d) - (h * h) < (r + size) * (r + size);//²¦¤ó©w²z
		}
		else {
			Vector3 d = (center - p2);
			return d.dot(d) < (r + size) * (r + size);
		}
	}
	bool collide(Vector3 center, float size) {
		return collide(center, size, *(new bool));
	}
	bool between(Vector3 center) {
		Vector3 vec = (p2 - p1);
		float ll = vec.dot(vec);
		float dot = vec.normalize().dot(center - p1);//vec normalized
		return dot>=0 && dot * dot <= ll;

	}

};

class Voxel {
public:
	Voxel() {};
	Voxel(Hash hash, Vector3 pos, int i) {
		this->hash = hash;
		this->pos = pos;
		this->i = i;
	}
	Hash hash;
	Hash hash2;
	Vector3 linkval;
	Vector3 pos;
	int i, ix, iy, iz;
	int belong;
	bool exist = true;
	bool immo = false;
	std::set<int> touchid;
	void setXYZ(int ix, int iy, int iz) {
		this->ix = ix;
		this->iy = iy;
		this->iz = iz;
	}
	void getXYZ(int& ix, int& iy, int& iz) {
		ix = this->ix;
		iy = this->iy;
		iz = this->iz;
	}
}; 
#endif // DATASTRUCTURE_H
