#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include "pch.h"
#include "mesh.h"

class Worth {
public:
	std::vector<float> val;
	Worth() {}
	Worth(std::vector<float> val_) {
		val = std::vector<float>(val_.size());
		for (int i = 0; i < val_.size(); i++) {
			val[i] = val_[i];
		}
	}
	bool operator < (const Worth& tar) const
	{
		for (int i = 0; i < val.size(); i++) {
			if (val[i] != tar.val[i])return val[i] < tar.val[i];
		}
		return false;
	}
	bool operator == (const Worth& tar) const
	{
		for (int i = 0; i < val.size(); i++) {
			if (val[i] != tar.val[i])return false;
		}
		return true;
	}
	bool operator > (const Worth& tar) const
	{
		for (int i = 0; i < val.size(); i++) {
			if (val[i] != tar.val[i])return val[i] > tar.val[i];
		}
		return false;
	}
};

class GroupLink {
public:
	GroupLink(int ida, int idb, Worth worth) {
		if (ida > idb) {
			this->ida = idb;
			this->idb = ida;
		}
		else {
			this->ida = ida;
			this->idb = idb;
		}
		this->worth = worth;
	}
	GroupLink(int ida, int idb) {
		if (ida > idb) {
			this->ida = idb;
			this->idb = ida;
		}
		else {
			this->ida = ida;
			this->idb = idb;
		}
	}
	int ida, idb;
	//float worth;
	Worth worth;
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
	TouchInfo(int e, Vector3 dir) {
		this->e = e;
		this->dir = dir;
	}
	TouchInfo(int e, Vector3 va, Vector3 vb, Vector3 dir, int knifeid) {
		this->e = e;
		this->va = va;
		this->vb = vb;
		this->dir = dir;
		this->knifeid = knifeid;
	}
	int e, knifeid;
	Vector3 va, vb, dir;
	bool operator < (const TouchInfo& tar) const
	{
		return e < tar.e;
	}
	bool operator == (const TouchInfo& tar) const
	{
		return e == tar.e;
	}
	bool operator > (const TouchInfo& tar) const
	{
		return e > tar.e;
	}
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
	Hash assignAll(int hashlen, bool b) {
		hash = std::vector<unsigned int>(hashlen / partlen + 1);
		if(!b)for (int i = 0; i < hash.size(); i++)hash[i] = 0x0;
		else  for (int i = 0; i < hash.size(); i++)hash[i] = 0xffffffff;
		return *this;
	}
	void assign(int i, bool b) {
		unsigned int mask = 0x1 << (i%partlen);
		if(b)hash[i / partlen] = hash[i / partlen] | mask;
		else hash[i / partlen] = hash[i / partlen] & (~mask);
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
	Hash operator|(const Hash& tar) {
		Hash newhash = *this;
		for (int i = 0; i < hash.size(); i++) {
			newhash.hash[i] = newhash.hash[i] | tar.hash[i];
		}
		return newhash;
	}
	Hash operator&(const Hash& tar) {
		Hash newhash = *this;
		for (int i = 0; i < hash.size(); i++) {
			newhash.hash[i] = newhash.hash[i] & tar.hash[i];
		}
		return newhash;
	}
	Hash operator^(const Hash& tar) {
		Hash newhash = *this;
		for (int i = 0; i < hash.size(); i++) {
			newhash.hash[i] = newhash.hash[i] ^ tar.hash[i];
		}
		return newhash;
	}
	int activeNum() {
		int cnt = 0;
		int lasthashlen = hashlen;
		for (int i = 0; i < hash.size(); i++, lasthashlen -= partlen)
		{
			unsigned int hashiter = hash[i];
			int curhashlen = lasthashlen < partlen ? lasthashlen : partlen;
			for (int j = 0; j < curhashlen; j++, hashiter /= 2) {
				if (hashiter % 2 == 0x1)cnt++;
			}
		}
		return cnt;
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

class Cube {
public:
	Vector3 ld, ru;
	Vector3 size;
	Vector3 cent;
	void apply(Vector3 ld, Vector3 ru) {
		this->ld = ld;
		this->ru = ru;
		size = (ru - ld);
		cent = (ru + ld) / 2;
	}
	Cube() {};
	Cube(Vector3 ld, Vector3 ru) {
		apply(ld, ru);
	}
	void mod(float val, int state) {
		if (state == 0)ru.x = std::min(ru.x, val);
		if (state == 1)ld.x = std::max(ld.x, val);
		if (state == 2)ru.y = std::min(ru.y, val);
		if (state == 3)ld.y = std::max(ld.y, val);
		if (state == 4)ru.z = std::min(ru.z, val);
		if (state == 5)ld.z = std::max(ld.z, val);
	}
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
	int it;//supervoxel
	std::vector<TouchInfo> touchinfos;
	std::vector<TouchInfo> boundinfos;
	std::set<int> boundids;
	std::vector<Vector3> voxels;
	std::vector<int> voxelsi;
	std::set<int> neighbor;//piece id
	std::map<int, std::vector<int>> contactmap;
	std::vector<int> upvoxels;
	std::vector<int> downvoxels;
	Vector3 max, min;
	int belong = -1;
	bool removed = false;
	bool iscube = false;
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
	std::vector<Vector3> avalist;
	std::vector<int> avalevel;
	std::set<int> avaset;
	std::map<int, std::vector<int>> contactmap;
	std::vector<std::map<int, std::vector<int>>> idxtcoldis = std::vector<std::map<int, std::vector<int>>>(6);
	std::vector<int> coldis;
	std::vector<int> borde;//voxels
	std::vector<std::map<int, int>> ztmap;
	std::map<int, std::set<int>> supervoxels;
	std::set<int> idxts;
	//bool isBound[6] = { false ,false ,false ,false ,false ,false };
	int avanum = 0;
	int bordenum = 0;
	//const float plus = 90.0f;
	//const float thre = 60.0f;
	const float plus = 45.0f;
	const float thre = 30.0f;
	const bool simpmode = true;// : utility.colmode
	bool trigger = false;
	bool removed = false;
	Vector3 maxDir;
	int maxDirDis;

	int id;
	std::set<int> neighbor;//id
	std::vector<int> upvoxels;
	std::vector<int> downvoxels;
	std::vector<std::set<int>> upplane;
	std::vector<std::set<int>> downplane;
	Vector3 max = Vector3(FLT_MIN, FLT_MIN, FLT_MIN), min = Vector3(FLT_MAX, FLT_MAX, FLT_MAX);
	bool zcav = false;
	float volume = 0;
	void initAva();
	void initAvaLev(int maxcol);
	Mesh getMesh(std::vector<Piece> & pieces);
	void operator=(const Group& group) {//useless?
		pieces = group.pieces;
		id = group.id;
		ava = group.ava;
		avalist = group.avalist;
		avanum = group.avanum;
		avaset = group.avaset;
		volume = group.volume;
		borde = group.borde;
		upplane = group.upplane;
		downplane = group.downplane;
		//contactmap = group.contactmap;
		idxtcoldis = group.idxtcoldis;
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
	float distance(Vector3 center) {
		Vector3 vec = (p2 - p1);
		float ll = vec.dot(vec);
		float dot = vec.normalize().dot(center - p1);//normalize
		float distance = 0;
		if (dot < 0) {
			Vector3 d = (center - p1);
			distance = d.length();
		}
		else if (dot*dot < ll) {
			Vector3 d = (center - p1);
			float h = vec.dot(d);
			distance = (d - h * vec).length();
		}
		else {
			Vector3 d = (center - p2);
			distance = d.length();
		}
		return distance > r ? distance - r : 0;
	}
	float distanceToSk(Vector3 center) {
		Vector3 vec = (p2 - p1);
		float ll = vec.dot(vec);
		float dot = vec.normalize().dot(center - p1);//normalize
		float distance = 0;
		if (dot < 0) {
			Vector3 d = (center - p1);
			distance = d.length();
		}
		else if (dot*dot < ll) {
			Vector3 d = (center - p1);
			float h = vec.dot(d);
			distance = (d - h * vec).length();
		}
		else {
			Vector3 d = (center - p2);
			distance = d.length();
		}
		return distance;
	}
	Vector3 proj(Vector3 center) {
		Vector3 vec = (p2 - p1);
		float ll = vec.dot(vec);
		float dot = vec.normalize().dot(center - p1);//normalize
		float distance = 0;
		if (dot < 0) {
			return p1;
		}
		else if (dot*dot < ll) {
			Vector3 d = (center - p1);
			float h = vec.dot(d);
			return h*vec + p1;
		}
		else {
			return p2;
		}
	}
	float rate(Vector3 pos){
		Vector3 proj = this->proj(pos);
		Vector3 vec = proj - p1;
		return vec.length() / (p2 - p1).length();
	}
	bool between(Vector3 center) {
		Vector3 vec = (p2 - p1);
		float ll = vec.dot(vec);
		float dot = vec.normalize().dot(center - p1);//vec normalized
		return dot>=0 && dot * dot <= ll;

	}
};

class Tag {
public:
	Tag() {}
	Tag(int e, int side, int lev){
		this->e = e;
		this->side = side;
		this->lev = lev;
	}
	int e;
	bool side;
	int lev;
	int tlev = 0;
	bool operator < (const Tag& tar) const
	{
		if (e == tar.e) {
			if (side == tar.side) {
				return lev < tar.lev;
			}
			return !side;
		}
		return e < tar.e;
	}
	bool operator == (const Tag& tar) const
	{
		return e == tar.e && side == tar.side && lev == tar.lev;
	}
	bool operator > (const Tag& tar) const
	{
		if (e == tar.e) {
			if (side == tar.side) {
				return lev > tar.lev;
			}
			return side;
		}
		return e > tar.e;
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
	Hash hash;//plane side (influence by hash2)
	Hash hash2;//related edge
	Hash hash3;//far collide
	Hash hash4;//manhatan 
	Tag tag;
	Vector3 linkval;
	float dis;
	Vector3 pos;
	int i, ix, iy, iz;
	int it, itx, ity, itz;
	int belong;
	bool exist = false;
	bool removed = false;
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
	bool tequal(const Voxel& tar) {
		return itx == tar.itx&&ity == tar.ity&&itz == tar.itz;
	}
	void setTXYZ(int itx, int ity, int itz) {
		this->itx = itx;
		this->ity = ity;
		this->itz = itz;
	}
	void getTXYZ(int& itx, int& ity, int& itz) {
		itx = this->itx;
		ity = this->ity;
		itz = this->itz;
	}
}; 


#endif // DATASTRUCTURE_H
