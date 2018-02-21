#include "utility.h"

//*******************************************
void Utility::voxelBfs() { //將所有touchedge先放入，並一組一組分別長
	std::set<int> vs;
	for (int i = 0; i < voxelsize; i++) vs.insert(i);
	for (int p = 0;; p++) {
		if(!voxelBfs(p, vs))break;
	}
	std::vector<std::vector<int>> queues;
	std::vector<int> flags;
	std::vector<std::set<int>> touchedges;
	for (int n = 0; n < pieces.size(); n++) {
		flags.push_back(0);
		queues.push_back(std::vector<int>());
		touchedges.push_back(std::set<int>());
		for (int i = 0; i < pieces[n].voxelsi.size(); i++) {
			int idx = pieces[n].voxelsi[i];
			queues[n].push_back(idx);
			std::vector<bool> bits = voxels[idx].hash.getBits();
			for (int j = 0; j < bits.size(); j += 3) {
				if (bits[j]) {
					touchedges[n].insert(j / 3);
				}
			}
		}
	}
	bool cont = true;
	while (cont) {
		cont = false;
		for (int n = 0; n < pieces.size(); n++) {
			voxelBfsOnce(n, queues[n], flags[n], touchedges[n], 0);
			if (flags[n] < queues[n].size())cont = true;
		}
	}
}

bool Utility::voxelBfs(int p, std::set<int> &vs) { //從未選的裡面找一個開始長
	std::vector<int> queue;
	int flag = 0;
	Piece piece;
	for (auto i : vs) {
		if (voxels[i].touchid.size() > 0 && !voxels[i].exist) {
			queue.push_back(i);
			voxels[i].exist = true;
			vs.erase(i);
			piece.voxelsi.push_back(i);
			voxels[i].belong = p;
			break;
		}
	}
	if (queue.size() == 0)return false;
	while (flag < queue.size()) {
		int qn = queue.size();
		for (; flag < qn; flag++) {
			int cur = queue[flag];

			int ix, iy, iz;
			voxels[cur].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			for (int i = 0; i < 6; i++) {
				if (nei[i] == -1)continue;
				int tar = nei[i];
				if(voxels[tar].exist)continue;
				int agreenum = 0;
				int disaggreenum = 0;
				for (int j = 0; j < topo.knifes.size();j++) {//only touch
					if (voxels[tar].hash.getBit(j * 3)) {
						if (voxels[cur].hash.getBit(j * 3 + 2) == voxels[tar].hash.getBit(j * 3 + 2))agreenum++;
						if (voxels[cur].hash.getBit(j * 3 + 2) != voxels[tar].hash.getBit(j * 3 + 2))disaggreenum++;
					}
				}
				if (agreenum > 0 && disaggreenum==0) {
					voxels[tar].exist = true;
					queue.push_back(tar);
					vs.erase(tar);
					piece.voxelsi.push_back(tar);
					voxels[tar].belong = p;
				}
			}
		}
	}
	if(piece.voxelsi.size()>8)pieces.push_back(piece); //purne small
	else {
		for (auto i : piece.voxelsi) {
			voxels[i].belong = -1;
		}
	}
	return true;
}

void Utility::voxelBfs(int p) {//從一個piece長到極限
	std::vector<int> queue;
	std::set<int> touchedge;
	int flag = 0;
	for (int i = 0; i < pieces[p].voxelsi.size(); i++) {
		int idx = pieces[p].voxelsi[i];
		queue.push_back(idx);
		std::vector<bool> bits = voxels[idx].hash.getBits();
		for (int j = 0; j < bits.size(); j += 3) {
			if (bits[j]) {
				touchedge.insert(j / 3);
			}
		}
	}
	while (flag < queue.size()) {
		int qn = queue.size();
		for (; flag < qn; flag++) {
			int cur = queue[flag];
			int ix, iy, iz;
			voxels[cur].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			for (int i = 0; i < 6; i++) {
				if (nei[i] == -1)continue;
				int tar = nei[i];
				if(voxels[tar].exist)continue;
				int agreenum = 0;
				int disaggreenum = 0;
				for (int j : touchedge) {
					if (voxels[cur].hash.getBit(j * 3 + 2) != voxels[tar].hash.getBit(j * 3 + 2)) disaggreenum++;
					else agreenum++;
				}
				if (agreenum > 0 && disaggreenum==0) {
					voxels[tar].exist = true;
					queue.push_back(tar);
					pieces[p].voxelsi.push_back(tar);
					voxels[tar].belong = p;
				}
			}
		}
	}
}

void Utility::voxelBfsOnce(int p, std::vector<int>& queue, int& flag, std::set<int>& touchedge, int mode) {//從一組開始長一次
	int qn = queue.size();
	for (; flag < qn; flag++) {
		int cur = queue[flag];
		int ix, iy, iz;
		voxels[cur].getXYZ(ix, iy, iz);
		int nei[6] = { -1, -1, -1, -1, -1, -1 };
		if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
		if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
		if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
		if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
		if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
		if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
		for (int i = 0; i < 6; i++) {
			if (nei[i] == -1)continue;
			int tar = nei[i];
			if (voxels[tar].exist)continue;
			int agreenum = 1;
			int disaggreenum = 0;
			for (int j : touchedge) {
				if (voxels[cur].hash.getBit(j * 3 + 2) != voxels[tar].hash.getBit(j * 3 + 2)) disaggreenum++;
				else agreenum++;
				/*
				if (voxels[cur].hash.getBit(j * 3 + 1)) {
					if (voxels[cur].hash.getBit(j * 3 + 2) != voxels[tar].hash.getBit(j * 3 + 2)) disaggreenum++;
					else agreenum++;
				}
				*/
			}
			if (agreenum > 0 && disaggreenum == 0) {
				voxels[tar].exist = true;
				queue.push_back(tar);
				pieces[p].voxelsi.push_back(tar);
				voxels[tar].belong = p;
			}
		}
	}
}

//*******************************************

int Utility::voxelId(int i, int j, int k, int mode) {
	int idx;
	if (mode == 0)idx = i * ny * nz + j * nz + k;
	if (mode == 1)idx = (nx - i - 1) * ny * nz + j * nz + k;
	if (mode == 2)idx = j * ny * nz + i * nz + k;
	if (mode == 3)idx = j * ny * nz + (ny - i - 1) * nz + k;
	if (mode == 4)idx = j * ny * nz + k * nz + i;
	if (mode == 5)idx = j * ny * nz + k * nz + (nz - i - 1);
	return idx;
}

int Utility::voxelTId(int i, int j, int k, int mode) {
	int idxt;
	if (mode == 0)idxt = i * nty * ntz + j * ntz + k;
	if (mode == 1)idxt = (ntx - i - 1) * nty * ntz + j * ntz + k;
	if (mode == 2)idxt = j * nty * ntz + i * ntz + k;
	if (mode == 3)idxt = j * nty * ntz + (nty - i - 1) * ntz + k;
	if (mode == 4)idxt = j * nty * ntz + k * ntz + i;
	if (mode == 5)idxt = j * nty * ntz + k * ntz + (ntz - i - 1);
	return idxt;
}

void Utility::voxelTId_(int & i, int & j, int & k, int mode, int idxt) {
	int itx = idxt / (nty * ntz);
	int ity = idxt / ntz % nty;
	int itz = idxt % ntz;
	if (mode == 0) { i = itx; j = ity; k = itz; }
	if (mode == 1) { i = ntx - itx - 1; j = ity; k = itz; }
	if (mode == 2) { j = itx; i = ity; k = itz; }
	if (mode == 3) { j = itx; i = nty - ity - 1; k = itz; }
	if (mode == 4) { j = itx; k = ity; i = itz; }
	if (mode == 5) { j = itx; k = ity; i = ntz - itz - 1; }
}

void Utility::genVoxelSeen() {//將moveable direction塞入hash，有計算linkval
	tic();
	//*****************************
	std::vector<bool> bits = std::vector<bool>(topo.knifes.size() * 3);
	nx = (ru.x - ld.x) / l + 1;
	ny = (ru.y - ld.y) / l + 1;
	nz = (ru.z - ld.z) / l + 1;
	voxels = new Voxel[voxelsize];
	//*******************************************
	#pragma omp parallel for
	for (int i = 0; i < voxelsize; i++) {
		int ix = i / (ny * nz);
		int iy = i / nz % ny;
		int iz = i % nz;
		Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
		int idx = ix * ny * nz + iy * nz + iz;
		voxels[idx] = Voxel(Hash(), pos, idx);
		voxels[idx].linkval = Vector3(0, 0, 0);
		bool immo = false;
		for (int j = 0; j < topo.knifes.size(); j++) {
			bool between;
			if (topo.caps[j].collide(pos, l / 2, between))immo = true;
			if (topo.caps[j].collide(pos, 1.5 *l, between))voxels[idx].touchid.insert(j);
		}
		voxels[idx].setXYZ(ix, iy, iz);
		voxels[idx].exist = false;
		voxels[idx].immo = immo;
	}
	//*******************************************
	#pragma omp parallel for
	for (int i = 0; i < voxelsize; i++) {
		int idx = i;
		int ix = voxels[idx].ix;
		int iy = voxels[idx].iy;
		int iz = voxels[idx].iz;
		Vector3 pos = voxels[idx].pos;
		voxels[idx].linkval = Vector3(FLT_MAX, FLT_MAX, FLT_MAX);
		int neiidxs[3] = { ix + 1 < nx ? voxelId(ix + 1, iy, iz, 0) : -1, iy + 1 < ny ? voxelId(ix, iy + 1, iz, 0) : -1, iz + 1 < nz ? voxelId(ix, iy, iz + 1, 0) : -1 };
		float maxlen = 100;
		for (int j = 0; j < topo.knifes.size(); j++) {
			float neivals[3];
			for (int k = 0; k < 3; k++) {
				if (neiidxs[k] != -1) {
					Vector3 tarpos = voxels[neiidxs[k]].pos;
					if (topo.knifes[j].distanceToPoint(pos)*topo.knifes[j].distanceToPoint(tarpos) < 0) {
						if (voxels[idx].touchid.count(j) > 0)neivals[k] = 0;
						else neivals[k] = (topo.knifes[j].center - (pos + tarpos) / 2).length();//暫時保持這樣，未來使用與cap距離
					}else neivals[k] = maxlen;//暫時保持這樣，未來需計算最大值
				}else neivals[k] = 0;
			}
			voxels[idx].linkval.x = std::min(voxels[idx].linkval.x, neivals[0]);
			voxels[idx].linkval.y = std::min(voxels[idx].linkval.y, neivals[1]);
			voxels[idx].linkval.z = std::min(voxels[idx].linkval.z, neivals[2]);
			
		}
	}
	//*******************************************
	for (int i = 0; i < 6; i++) {
		calVoxelSeen(i);
	}
	//*******************************************
	std::set<Hash> hs;
	for (int i = 0; i < voxelsize; i++)hs.insert(voxels[i].hash);
	std::cout << "Label number : " << hs.size() << std::endl;
	//*******************************************
	toc();
}

void Utility::calVoxelSeen(int mode) {//計算一個方向的塞入
	int n, w1, w2;
	Vector3 dir;
	if (mode/2==0) { n = nx; w1 = ny; w2 = nz; dir = Vector3(-1, 0, 0); }
	if (mode/2==1) { n = ny; w1 = nx; w2 = nz; dir = Vector3(0, -1, 0); }
	if (mode/2==2) { n = nz; w1 = nx; w2 = ny; dir = Vector3(0, 0, -1); }
	if (mode%2==1) dir *= -1;
	std::vector<int> pixel(w1*w2); for (int j = 0 ; j < pixel.size(); j++) pixel[j] = 0;
	int cnt = 0;

	for (int i = 0; i < n; i++)  {
		for (int j = 0; j < w1; j++) {
			for (int k = 0; k < w2; k++) {
				int idx = voxelId(i, j, k, mode);
				if (voxels[idx].immo) pixel[j * w2 + k] = 1;
				bool crit = false;
				for (auto t : voxels[idx].touchid) {
					Vector3 sampsir = topo.knifes[t].normal.normalize();
					if (topo.knifes[t].distanceToPoint(voxels[idx].pos) < 0) sampsir *= -1;
					//if (state==3 && topo.knifes[t].distanceToPoint(voxels[idx].pos) < 0) std::cout << dir.dot(sampsir) << " " << cos(105 * M_PI / 180) << std::endl;
					if (dir.dot(sampsir) < cos(75 * M_PI / 180)) crit = true;
				}
				/**///下方全部刪除，因為有反面的檢查，所以效果一樣?
				if (crit) pixel[j * w2 + k] = 1;
				/**/
				if (pixel[j * w2 + k] == 0 && !crit) voxels[idx].hash.addHash(true);
				//if (pixel[j * w2 + k] == 0) voxels[idx].hash.addHash(true);
				else voxels[idx].hash.addHash(false);

				//if (voxels[idx].immo)cnt++;
				if (pixel[j * w2 + k] == 0)cnt++;
			}
		}
	}
	std::cout << "Mode " << mode << " : " << cnt << std::endl;
}

void Utility::voxelDirSearch(int mode) {
	Piece piece;
	int p = pieces.size();
	int n, w1, w2;
	//*******************************************
	if (mode / 2 == 0) { n = nx; w1 = ny; w2 = nz; }
	if (mode / 2 == 1) { n = ny; w1 = nx; w2 = nz; }
	if (mode / 2 == 2) { n = nz; w1 = nx; w2 = ny; };

	//*******************************************
	std::vector<int> pixel(w1*w2); for (int j = 0; j < w1; j++) for (int k = 0; k < w2; k++) pixel[j*w2 + k] = 0;
	bool nfirst = false;
	for (int i = 0; i < n; i++) {
		nfirst = true;
		int anything = 0;
		for (int j = 0; j < w1; j++) {
			for (int k = 0; k < w2; k++) {
				int idx = voxelId(i, j, k, mode);
				if (voxels[idx].immo) pixel[j*w2 + k] = 1;
				if (pixel[j*w2 + k] == 0 && !voxels[idx].exist) {
					piece.voxelsi.push_back(idx);
					voxels[idx].exist = true;
					voxels[idx].belong = p;
					anything++;
					nfirst = true;
				}
				/*
				int idx_prev = voxelId(i-1, j, k, mode);
				if (!voxels[idx].immo && !voxels[idx].exist) {
				if (!nfirst) {
				piece.voxelsi.push_back(idx);
				voxels[idx].exist = true;
				voxels[idx].belong = p;
				anything++;
				}
				else if (voxels[idx_prev].exist && voxels[idx_prev].belong == p) {
				piece.voxelsi.push_back(idx);
				voxels[idx].exist = true;
				voxels[idx].belong = p;
				anything++;
				}
				}
				*/
			}
		}
		if (anything == 0 && !nfirst)break;
	}
	pieces.push_back(piece);

	//*******************************************
	Group group = Group();
	if (mode == 0) { group.ava.push_back(Vector3(-1, 0, 0)); }
	if (mode == 2) { group.ava.push_back(Vector3(0, -1, 0)); }
	if (mode == 4) { group.ava.push_back(Vector3(0, 0, -1)); }
	if (mode == 1) { group.ava.push_back(Vector3(1, 0, 0)); }
	if (mode == 3) { group.ava.push_back(Vector3(0, 1, 0)); }
	if (mode == 5) { group.ava.push_back(Vector3(0, 0, 1)); }
	group.avanum = 1;
	groups.push_back(group);
	//*******************************************
}

void Utility::voxelDirSeen() {
	voxelDirSeen(0);
	voxelDirSeen(1);
	voxelDirSeen(2);
	voxelDirSeen(3);
	voxelDirSeen(4);
	voxelDirSeen(5);
}

/*
void Utility::voxelDirSeen(int mode) {
	Piece piece;
	int p = pieces.size();
	int n, w1, w2;
	//*******************************************
	if (mode / 2 == 0) { n = nx; w1 = ny; w2 = nz; }
	if (mode / 2 == 1) { n = ny; w1 = nx; w2 = nz; }
	if (mode / 2 == 2) { n = nz; w1 = nx; w2 = ny; }
	//*******************************************
	std::vector<int> pixel(w1*w2);
	for (int j = 0; j < pixel.size(); j++) pixel[j] = 0;
	bool nfirst = false;
	float jj = ((float)w1 / 2) / ((float)n / 2);
	float kk = ((float)w2 / 2) / ((float)n / 2);
	for (int i = 0; i < (n+1)/2; i++, nfirst = true) {
		nfirst = true;
		int anything = 0;
		for (int j = jj*i; j < w1 - jj * i; j++) {
			for (int k = kk*i; k < w2 - kk * i; k++) {
				int idx = voxelId(i, j, k, mode);
				if (voxels[idx].immo) pixel[j*w2 + k] = 1;
				if (pixel[j*w2 + k] == 0 && !voxels[idx].exist) {
					piece.voxelsi.push_back(idx);
					voxels[idx].exist = true;
					voxels[idx].belong = p;
					anything++;
				}
			}
		}
		if (anything == 0)break;
	}
	pieces.push_back(piece);
	//*******************************************
	Group group = Group();
	if (mode == 0) { group.ava.push_back(Vector3(-1, 0, 0)); }
	if (mode == 2) { group.ava.push_back(Vector3(0, -1, 0)); }
	if (mode == 4) { group.ava.push_back(Vector3(0, 0, -1)); }
	if (mode == 1) { group.ava.push_back(Vector3(1, 0, 0)); }
	if (mode == 3) { group.ava.push_back(Vector3(0, 1, 0)); }
	if (mode == 5) { group.ava.push_back(Vector3(0, 0, 1)); }
	group.avanum = 1;
	groups.push_back(group);
	//*******************************************
}
*/

void Utility::voxelDirSeen(int mode) {
	Piece piece;
	int p = pieces.size();
	int n, w1, w2;
	//*******************************************
	if (mode / 2 == 0) { n = nx; w1 = ny; w2 = nz; }
	if (mode / 2 == 1) { n = ny; w1 = nx; w2 = nz; }
	if (mode / 2 == 2) { n = nz; w1 = nx; w2 = ny; }
	//*******************************************
	std::vector<int> pixel(w1*w2);
	for (int j = 0; j < pixel.size(); j++) pixel[j] = 0;
	for (int i = 0; i < n; i++) {
		#pragma omp parallel for
		for (int j = 0; j < w1; j++) {
			for (int k = 0; k < w2; k++) {
				int idx = voxelId(i, j, k, mode);
				if (voxels[idx].dis < l) pixel[j * w2 + k] = 1;
				if (pixel[j*w2 + k] == 0) {
					voxels[idx].hash3.addHash(true);
				}
				else {
					voxels[idx].hash3.addHash(false);
				}
			}
		}
	}
	//*******************************************
}

void Utility::voxelBfsSeen() {
	std::set<int> last;
	for (int i = 0; i < voxelsize; i++)last.insert(i);

	while (voxelBfsSeenOnce(last, 4));
	while (voxelBfsSeenOnce(last, 5));
	//voxelCollectSeen(0);
	//voxelCollectSeen(1);
	//voxelCollectSeen(2);
	//voxelCollectSeen(3);
	//while (voxelBfsSeenOnce(last, 0));
	//while (voxelBfsSeenOnce(last, 1));
	//while (voxelBfsSeenOnce(last, 2));
	//while (voxelBfsSeenOnce(last, 3));
}

bool Utility::voxelBfsSeenOnce(std::set<int> & last, int mode) {

	std::vector<int> queue;
	int flag = 0;
	Piece piece;
	int p = pieces.size();
	float min = FLT_MAX;int mini;
	for (auto i : last) {
		if (voxels[i].immo || voxels[i].exist)continue;
		if (voxels[i].touchid.size() == 0)continue;//從touch出發
		if (!voxels[i].hash.getBit(mode))continue;
		int ix, iy, iz;
		voxels[i].getXYZ(ix, iy, iz);
		if (mode == 0)if (min > ix) { min = ix; mini = i; };
		if (mode == 1)if (min > -ix) { min = -ix; mini = i; };
		if (mode == 2)if (min > iy) { min = iy;  mini = i;};
		if (mode == 3)if (min > -iy) { min = -iy;  mini = i;};
		if (mode == 4)if (min > iz) { min = iz;  mini = i;};
		if (mode == 5)if (min > -iz) { min = -iz;  mini = i;};
	}
	if (min == FLT_MAX )return false;
	queue.push_back(mini);
	piece.voxelsi.push_back(mini);
	voxels[mini].exist = true;
	voxels[mini].belong = p;
	while (flag < queue.size()) {
		int qn = queue.size();
		for (; flag < qn; flag++) {
			int cur = queue[flag];
			last.erase(cur);
			int ix, iy, iz;
			voxels[cur].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			for (int i = 0; i < 6; i++) {
				int tar = nei[i];
				if (tar == -1)continue;
				if (voxels[tar].exist)continue;
				if(voxels[tar].immo)continue;
				//if (voxels[tar].touchid.size() > 0)
				if (voxels[tar].hash.getBit(mode)) {
					queue.push_back(tar);
					piece.voxelsi.push_back(tar);
					voxels[tar].exist = true;
					voxels[tar].belong = p;
				}
				
			}
		}
	}
	//*******************************************
	
	//if (piece.voxelsi.size() < 100) {//暫時修改
	if (false) {
		for (int i = 0; i < piece.voxelsi.size(); i++) {
			voxels[piece.voxelsi[i]].exist = false;
		}
	}
	else {
		pieces.push_back(piece);
		genVoxelMesh(p, mode);
		Group group = Group();
		if (mode == 0) { group.ava.push_back(Vector3(-1, 0, 0)); }
		if (mode == 2) { group.ava.push_back(Vector3(0, -1, 0)); }
		if (mode == 4) { group.ava.push_back(Vector3(0, 0, -1)); }
		if (mode == 1) { group.ava.push_back(Vector3(1, 0, 0)); }
		if (mode == 3) { group.ava.push_back(Vector3(0, 1, 0)); }
		if (mode == 5) { group.ava.push_back(Vector3(0, 0, 1)); }
		group.avanum = group.ava.size();
		groups.push_back(group);
	}
	//*******************************************
}

//*******************************************

void Utility::voxelCollectSeen(int mode) {//group只用來記錄方向
	Piece piece;
	for (int i = 0; i < voxelsize; i++) {
		if (!voxels[i].exist && !voxels[i].immo && voxels[i].hash.getBit(mode)) {
			piece.voxelsi.push_back(i);
			voxels[i].exist = true;
		}
	}
	pieces.push_back(piece);
	//***************************
	Group group = Group();
	if (mode == 0) { group.ava.push_back(Vector3(-1, 0, 0)); }
	if (mode == 2) { group.ava.push_back(Vector3(0, -1, 0)); }
	if (mode == 4) { group.ava.push_back(Vector3(0, 0, -1)); }
	if (mode == 1) { group.ava.push_back(Vector3(1, 0, 0)); }
	if (mode == 3) { group.ava.push_back(Vector3(0, 1, 0)); }
	if (mode == 5) { group.ava.push_back(Vector3(0, 0, 1)); }
	group.avanum = group.ava.size();
	groups.push_back(group);
}

void Utility::collectLast() {
	Piece piece;
	for (int i = 0; i < voxelsize; i++) {
		if (!voxels[i].exist && !voxels[i].immo) {
			voxels[i].exist = true;
			piece.voxelsi.push_back(i);
			voxels[i].belong = pieces.size();
		}
	}
	if (piece.voxelsi.size() > 0) {
		pieces.push_back(piece);
		groups.push_back(Group());
	}
}

void Utility::purneVoxel(Group & group) {
	int g = group.id;
	for (int nn = 0; nn < group.pieces.size(); nn++) {
		int n = group.pieces[nn];
		std::vector<int> voxelsi_;
		for (int m = 0; m < pieces[n].voxelsi.size(); m++) {
			int i = pieces[n].voxelsi[m];
			if (!voxels[i].exist) continue;
			int ix, iy, iz;
			voxels[i].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1) * nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1) * nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			int samecnt = 0;
			for (int j = 0; j < 6; j++) {
				if (nei[j] == -1)continue;
				if (!voxels[nei[j]].exist)continue;
				if (voxels[nei[j]].belong != -1 && voxels[nei[j]].belong != g)continue;
				samecnt++;
			}
			if (samecnt < 6) {
				voxelsi_.push_back(i);
				
			}
			else voxels[i].belong == -1;
		}
		pieces[n].voxelsi.swap(voxelsi_);
	}
}

void Utility::genVoxelOutput() {
	for (int n = 0; n < pieces.size(); n++) {
		for (int m = 0; m < pieces[n].voxelsi.size(); m++) {
			int i = pieces[n].voxelsi[m];
			if (!voxels[i].exist) continue;
			int ix, iy, iz;
			voxels[i].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1) * nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1) * nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			int samecnt = 0;
			for (int j = 0; j < 6; j++) {
				if (nei[j] == -1)continue;
				if (!voxels[nei[j]].exist)continue;
				if (voxels[nei[j]].belong != -1 && voxels[nei[j]].belong != n)continue;
				samecnt++;
			}
			if (samecnt < 6) {
				pieces[n].voxels.push_back(voxels[i].pos);
			}
		}
	}
}

void Utility::genVoxelMesh(int p, int mode) {
	Piece piece = pieces[p];
	int n, w1, w2;
	//*******************************************
	if (mode/2==0 ) { n = nx; w1 = ny; w2 = nz; }
	if (mode/2==1) { n = ny; w1 = nx; w2 = nz; }
	if (mode/2==2) { n = nz; w1 = nx; w2 = ny; }
	//*******************************************
	std::vector<int> pixelmin(w1*w2);
	std::vector<int> pixelmax(w1*w2);
	std::vector<int> pixel(w1*w2);
	for (int j = 0; j < pixelmin.size(); j++) pixelmin[j] = 0;
	for (int j = 0; j < pixelmax.size(); j++) pixelmax[j] = n-1;
	for (int j = 0; j < pixel.size(); j++) pixel[j] = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < w1; j++) {
			for (int k = 0; k < w2; k++) {
				int idx = voxelId(i, j, k, mode);
				int idx_ = voxelId(n-i-1, j, k, mode);
				if (voxels[idx_].exist && voxels[idx_].belong == p) {
					pixelmin[j*w2 + k] = n - i - 1;
				}
				if (voxels[idx].exist && voxels[idx].belong == p) {
					pixelmax[j * w2 + k] = i;
					pixel[j * w2 + k] = 1;
				}
			}
		}
	}
	Mesh mesh1;
	for (int j = 0; j < w1; j++) {
		for (int k = 0; k < w2; k++) {
			int i = pixelmin[j * w2 + k];
			int idx = voxelId(i, j, k, mode);
			mesh1.pushVertice(voxels[idx].pos);
			if (j < w1 - 1 && k < w2 - 1) {
				int a = j * w2 + k;
				int b = j * w2 + (k + 1);
				int c = (j + 1) * w2 + (k + 1);
				int d = (j + 1) * w2 + k;
				/*
				if (state > 0)mesh1.pushIndice(a, b, c, d);
				else mesh1.pushIndice(a, d, c, b);
				*/
				if (pixel[a] == 1) {
					if (pixel[b] == 1 && pixel[c] == 1) {
						if(mode %2 == 0)mesh1.pushIndice(a, b, c);
						else mesh1.pushIndice(a, c, b);
					}
					if (pixel[c] == 1 && pixel[d] == 1) {
						if (mode % 2 == 0)mesh1.pushIndice(a, c, d);
						else mesh1.pushIndice(a, d, c);
					}
				}
			}
		}
	}
	Mesh mesh2;
	for (int j = 0; j < w1; j++) {
		for (int k = 0; k < w2; k++) {
			int i = pixelmax[j * w2 + k];
			int idx = voxelId(i, j, k, mode);
			mesh2.pushVertice(voxels[idx].pos);
			if (j < w1 - 1 && k < w2 - 1) {
				int a = j * w2 + k;
				int b = j * w2 + (k + 1);
				int c = (j + 1) * w2 + (k + 1);
				int d = (j + 1) * w2 + k;
				if (pixel[a] == 1) {
					if (pixel[b] == 1 && pixel[c] == 1) {
						if (mode % 2 == 0)mesh2.pushIndice(a, c, b);
						else mesh2.pushIndice(a, b, c);
					}
					if (pixel[c] == 1 && pixel[d] == 1) {
						if (mode % 2 == 0)mesh2.pushIndice(a, d, c);
						else mesh2.pushIndice(a, c, d);
					}
				}
			}
		}
	}
	mesh1.addMesh(mesh2);
	int prefix = w1 * w2;
	std::vector<int> pixeledge(w1 * w2);
	for (int j = 0; j < pixel.size(); j++) pixeledge[j] = 0;
	for (int j = 0; j < w1; j++) {
		for (int k = 0; k < w2; k++) {
			if (pixel[j * w2 + k]) {
				if (j - 1 < 0 || pixel[(j - 1) * w2 + k] == 0)pixeledge[j * w2 + k] = 1;
				if (j + 1 >= w1 || pixel[(j + 1) * w2 + k] == 0)pixeledge[j * w2 + k] = 1;
				if (k - 1 < 0 || pixel[j * w2 + (k - 1)] == 0)pixeledge[j * w2 + k] = 1;
				if (k + 1 >= w2 || pixel[j * w2 + (k + 1)] == 0)pixeledge[j * w2 + k] = 1;
			}
		}
	}
	for (int j = 0; j < w1; j++) {
		for (int k = 0; k < w2; k++) {
			int a = j * w2 + k;
			int b = j * w2 + (k + 1);
			int d = (j + 1) * w2 + k;
			if (pixeledge[a] == 1) {
				if (k + 1 < w2 && pixeledge[b] == 1) {
					if (mode % 2 == 0) {
						if (j + 1 >= w1 || pixel[(j + 1) * w2 + k] == 0)
							mesh1.pushIndice(a, b, b + prefix, a + prefix);
						if (j - 1 < 0 || pixel[(j - 1) * w2 + k] == 0)
							mesh1.pushIndice(a, a + prefix, b + prefix, b);
					}
					else {
						if (j + 1 >= w1 || pixel[(j + 1) * w2 + k] == 0)
							mesh1.pushIndice(a, a + prefix, b + prefix, b);
						if (j - 1 < 0 || pixel[(j - 1) * w2 + k] == 0)
							mesh1.pushIndice(a, b, b + prefix, a + prefix);
					}
					
				}
				if (j + 1 < w1 && pixeledge[d] == 1) {
					if (mode%2 == 0) {
						if (k + 1 >= w2 || pixel[j * w2 + (k + 1)] == 0)
							mesh1.pushIndice(a, a + prefix, d + prefix, d);
						if (k - 1 < 0 || pixel[j * w2 + (k - 1)] == 0)
							mesh1.pushIndice(a, d, d + prefix, a + prefix);
					}
					else {
						if (k + 1 >= w2 || pixel[j * w2 + (k + 1)] == 0)
							mesh1.pushIndice(a, d, d + prefix, a + prefix);
						if (k - 1 < 0 || pixel[j * w2 + (k - 1)] == 0)
							mesh1.pushIndice(a, a + prefix, d + prefix, d);
					}
				}
			}
		}
	}
	pieces[p].mesh = mesh1;
	pieces[p].mesh.purne();
}

//*******************************************

void Utility::genVoxel() {
	tic();
	std::vector<bool> bits = std::vector<bool>(topo.knifes.size() * 3);
	nx = (ru.x - ld.x) / l + 1;
	ny = (ru.y - ld.y) / l + 1;
	nz = (ru.z - ld.z) / l + 1;
	//voxels = std::vector<Voxel>(nx*ny*nz);
	voxels = new Voxel[voxelsize];
	#pragma omp parallel for
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
				int idx = ix * ny * nz + iy * nz + iz;
				voxels[idx] = Voxel(Hash(), pos, idx);
				std::vector<bool> bits = std::vector<bool>(topo.knifes.size()*3);
				bool immo = false;
				for (int i = 0; i < topo.knifes.size(); i++) {
					bool between;
					bits[i * 3] = topo.caps[i].collide(pos, l, between);//it is L
					bits[i * 3 + 1] = between;
					bits[i * 3 + 2] = (pos - topo.knifes[i].center).dot(topo.knifes[i].normal) >= 0;
					if (bits[i * 3]) voxels[idx].touchid.insert(i);
					if (topo.caps[i].collide(pos, l/2, between))immo = true;
				}
				voxels[idx].hash.assign(bits);
				voxels[idx].setXYZ(ix, iy, iz);
				voxels[idx].exist = false;
				voxels[idx].immo = immo;
			}
		}
	}
	toc();
	//**************************************
	//test
	/* 
	pieces = std::vector<Piece>(topo.edgenum);
	for (int i = 0; i < voxelsize; i++) {
		std::vector<bool> bits = voxels[i].hash.getBits();
		for (int j = 0; j < bits.size(); j+=3) {
			if (bits[j] && bits[j+2]) {
				pieces[j / 3].voxels.push_back(voxels[i].pos);
				pieces[j / 3].voxelsi.push_back(i);
				pieces[j / 3].volume++;
				break;
			}
		}
	}
	for (int i = 0; i < pieces.size(); i++) {
		std::cout << pieces[i].voxelsize << std::endl;
	}
	*/
}

void Utility::genVoxelByKnife_autotune() {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> rand(-5, 5);
	std::uniform_real_distribution<float> rand2(10, 15);
	float min = FLT_MAX;
	Vector3 minmove;
	for (int i = 0; i < 30; i++) {
		std::cout << "Attemp " << i << std::endl;
		if (i > 0) {
			move = Vector3(rand(generator), rand(generator), rand(generator));
			//stx = sty = stz = std::round(rand2(generator));
		}
		calBound();
		float val = genVoxelByKnife();
		if (val < min) {
			min = val;
			minmove = move;
			std::cout << move << " " << stx << " " << val << std::endl;
		}
	}
	move = minmove;
	calBound();
	genVoxelByKnife();
}

void Utility::genSuperVoxelBlock() {
	//voxels = std::vector<Voxel>(nx*ny*nz);//calBound has do
	//***************************************************************
	tic();
	std::vector<std::set<int>> tnear(ntx * nty * ntz);
	topo.tnear = std::vector<std::set<int>>(ntx * nty * ntz);
	float ltx = stx * l / 2, lty = sty * l / 2, ltz = stz * l / 2;
	//float lt = std::sqrtf(ltx*ltx + lty*lty + ltz*ltz);
	float lt = std::min(ltx, std::min(lty,ltz));//short mode
	float llt = sqrt(ltx * ltx + lty * lty + ltz * ltz);
	#pragma omp parallel for
	for (int itx = 0; itx < ntx; itx++) {
		for (int ity = 0; ity < nty; ity++) {
			for (int itz = 0; itz < ntz; itz++) {
				Vector3 pos = ld + Vector3(l * (stx - 1) / 2, l * (sty - 1) / 2, l * (stz - 1) / 2) + Vector3(l * stx * itx, l * sty * ity, l * stz * itz);
				int idxt = itx * nty * ntz + ity * ntz + itz;
				for (int i = 0; i < topo.knifes.size(); i++) {
					if (topo.caps[i].collide(pos, lt))tnear[idxt].insert(i);
					if (topo.caps[i].collide(pos, llt))topo.tnear[idxt].insert(i);
				}
			}
		}
	}
	printf("cal near - "); toc();
	//***************************************************************
	supervoxelblock = std::vector<std::vector<int>>(ntx * nty * ntz, std::vector<int>(6, 1));// : free ; 0 : block
	topo.edgeblock = std::vector<std::vector<int>>(topo.edgenum, std::vector<int>(6, 1));// : free ; 0 : block
	for (int mode = 0; mode < 6; mode++) {
		int n, w1, w2;
		if (mode / 2 == 0) { n = ntx; w1 = nty; w2 = ntz;}
		if (mode / 2 == 1) { n = nty; w1 = ntx; w2 = ntz;}
		if (mode / 2 == 2) { n = ntz; w1 = ntx; w2 = nty;}
		int mode_ = (mode / 2 * 2) + ((mode % 2 + 1) % 2);
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < w1; j++) {
				for (int k = 0; k < w2; k++) {
					int tidxp = voxelTId(i-1, j, k, mode);
					int tidx = voxelTId(i, j, k, mode);
					if (tnear[tidxp].size() == 0) {
						supervoxelblock[tidx][mode] = supervoxelblock[tidxp][mode];
					}
					else {
						supervoxelblock[tidx][mode] = 0;
					}
				}
			}
		}
		for (int ti = 0; ti < tnear.size(); ti++) {
			for (auto e : tnear[ti]) {
				topo.edgeblock[e][mode] = std::min(topo.edgeblock[e][mode], supervoxelblock[ti][mode]);
			}
		}
	}
	for (int e = 0; e < topo.edgeblock.size(); e++) {
		printf("Edge %d : ", e);
		for (int i = 0; i < 6; i++) {
			printf("%d ", topo.edgeblock[e][i]);
		}
		printf("\n");
	}
	
}

void Utility::preview() {
	previewVoxelByKnife(topo.knifes);
	previewPiece_voxel(topo.knifes);
}

void Utility::previewVoxelByKnife(std::vector<Plane>& knifes) {
	tbits2 = std::vector<std::vector<bool>>(ntx * nty * ntz, std::vector<bool>(knifes.size(), false));
	//*******************************************************************
	tic();
	float ltx = stx * l / 2, lty = sty * l / 2, ltz = stz * l / 2;
	float lt = std::sqrtf(ltx*ltx + lty*lty + ltz*ltz);
	lt = std::min(std::min(ltx, lty), ltz);
	#pragma omp parallel for
	for (int itx = 0; itx < ntx; itx++) {
		for (int ity = 0; ity < nty; ity++) {
			for (int itz = 0; itz < ntz; itz++) {
				Vector3 pos = ld + Vector3(l * (stx - 1) / 2, l * (sty - 1) / 2, l * (stz - 1) / 2) + Vector3(l * stx * itx, l * sty * ity, l * stz * itz);
				int idxt = itx * nty * ntz + ity * ntz + itz;
				for (int i = 0; i < knifes.size(); i++) {
					if (topo.caps[i].collide(pos, lt)) {
						tnear[idxt].push_back(i);
						tbits2[idxt][i] = true;
					}
				}
			}
		}
	}
	if (printdebug) { printf("cal near - "); toc(); }
	//*******************************************************************
	tic();
	#pragma omp parallel for
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
				int idx = ix * ny * nz + iy * nz + iz;
				int itx = ix / stx, ity = iy / sty, itz = iz / stz;
				int idxt = itx * nty * ntz + ity * ntz + itz;
				voxels[idx] = Voxel(Hash(), pos, idx);
				voxels[idx].setXYZ(ix, iy, iz);
				voxels[idx].setTXYZ(itx, ity, itz);
				voxels[idx].it = idxt;
				voxels[idx].hash.assignAll(knifes.size(), false);
			}
		}
	}
	if (printdebug) { printf("cal vox  - "); toc(); }
	//*******************************************************************
}

void Utility::previewPiece_voxel(std::vector<Plane>& knifes) {
	tic();
	//*******************************************************************
	#pragma omp parallel for
	for (int idx = 0; idx<voxelsize; idx++) {
		int idxt = voxels[idx].it;
		Vector3 pos = voxels[idx].pos;
		for (auto i : tnear[idxt]) {
			bool bit = knifes[i].distanceToPoint(pos) >= 0;
			voxels[idx].hash.assign(i, bit & tbits2[idxt][i]);
		}
	}
	if (printdebug) { printf("cal hash - "); toc(); }
	//*******************************************************************
	tic();
	std::vector<std::map<Hash, int>> hashmaps(ntx*nty*ntz);
	for(int idxt = 0; idxt < tnear.size(); idxt++){
		//int idxt = itx * nty * ntz + ity * ntz + itz;
		int itx = idxt / (nty*ntz);
		int ity = idxt / ntz % nty;
		int itz = idxt % ntz;
		//printf("%d %d %d : %d %d %d\n", itx, ity, itz, ntx, nty, ntz);
		for (int ix = itx * stx; ix < itx * stx + stx; ix++) for (int iy = ity * sty; iy < ity * sty + sty; iy++) for (int iz = itz * stz; iz < itz * stz + stz; iz++) {
			int idx = ix * ny * nz + iy * nz + iz;
			Hash hash = voxels[idx].hash;
			if (hashmaps[idxt].count(hash) == 0) {
				Piece piece;
				piece.id = pieces.size();
				hashmaps[idxt][hash] = piece.id;
				pieces.push_back(piece);
			}
			int p = hashmaps[idxt][hash];
			voxels[idx].belong = p;
			voxels[idx].exist = true;
			pieces[p].voxelsi.push_back(idx);
			pieces[p].volume++;
		}
	}
	if (printdebug) { printf("cal pn  - ");  toc();}
	tic();
	#pragma omp parallel for
	for (int p = 0; p < pieces.size(); p++) {
		Voxel vox = voxels[pieces[p].voxelsi[0]];
		int idxt = vox.it;
		for (auto e : tnear[idxt]) {
			pieces[p].touchinfos.push_back(TouchInfo(e, knifes[e].normal * (vox.hash.getBit(e) ? 1 : -1)));
		}
		pieces[p].id = p;
		pieces[p].it = vox.it;
	}
	if (printdebug) { printf("cal ti   - ");  toc(); printf("Generate Piece %d\n", pieces.size());}
	
	//*****************************************************************
	tic();
	#pragma omp parallel for
	for (int p = 0; p < pieces.size(); p++) {
		topo.calTouchBound(pieces[p], false);
		if(pieces[p].boundids.size()==topo.avalist.size())
			topo.calTouchBound(pieces[p], true);
	}
	if (printdebug) { printf("cal bound- "); toc(); }
}

float Utility::genVoxelByKnife() {
	//voxels = std::vector<Voxel>(nx*ny*nz);
	//***************************************************************
	tic();
	std::vector<float> ratefrom(topo.knifes.size());
	std::vector<float> rateto(topo.knifes.size());
	#pragma omp parallel for
	for (int i = 0; i < topo.knifes.size(); i++) {
		Vector3 v1 = topo.vertices[topo.edges[i].ia];
		Vector3 v2 = topo.vertices[topo.edges[i].ib];
		Vector3 vec = (v2 - v1).normalize();
		ratefrom[i] = topo.caps[i].rate(v1 + vec * topo.edges[i].fixa);
		rateto[i] = topo.caps[i].rate(v2 - vec * topo.edges[i].fixb);
		ratefrom[i] = ratefrom[i] > 0.4 ? 0.4 : ratefrom[i];
		rateto[i] = rateto[i] < 0.6 ? 0.6 : rateto[i];
		//std::cout << ratefrom[i] << " ";
		//std::cout << rateto[i] << std::endl;
	}
	if (printdebug) { printf("cal rate - "); toc(); }
	//*******************************************************************
	tic();
	std::vector<std::vector<int>> tnear(ntx * nty * ntz);
	float ltx = stx * l / 2, lty = sty * l / 2, ltz = stz * l / 2;
	float lt = std::sqrtf(ltx*ltx + lty*lty + ltz*ltz);
	#pragma omp parallel for
	for (int itx = 0; itx < ntx; itx++) {
		for (int ity = 0; ity < nty; ity++) {
			for (int itz = 0; itz < ntz; itz++) {
				Vector3 pos = ld + Vector3(l * (stx - 1)  / 2, l * (sty - 1) / 2, l * (stz - 1) / 2) + Vector3(l * stx * itx, l * sty * ity, l * stz * itz);
				int idxt = itx * nty * ntz + ity * ntz + itz;
				for (int i = 0; i < topo.knifes.size(); i++) {
					if(topo.caps[i].collide(pos, lt))tnear[idxt].push_back(i);
				}
			}
		}
	}
	if (printdebug) {printf("cal near - "); toc();}
	//***************************************************************
	tic();
	#pragma omp parallel for
	for (int ix = 0; ix < nx; ix++) {
		//printf("xLayer %d\n", ix);
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
				int idx = ix * ny * nz + iy * nz + iz;
				int itx = ix / stx, ity = iy / sty, itz = iz / stz;
				int idxt = itx * nty * ntz + ity * ntz + itz;
				std::vector<bool> bits = std::vector<bool>(topo.knifes.size());
				float mincentdis = FLT_MAX;
				int mincent = -1;
				float dis = FLT_MAX;
				voxels[idx] = Voxel(Hash(), pos, idx);
				for(auto i : tnear[idxt]){
				//for (int i = 0; i < topo.knifes.size(); i++) {
					bits[i] = topo.knifes[i].distanceToPoint(pos) >= 0;
					float curdis = topo.caps[i].distance(pos);
					float curdisToSk = topo.caps[i].distanceToSk(pos);
					//bool between, collide = topo.caps[i].collide(pos, l, between);
					float rate = topo.caps[i].rate(pos);
					if (curdisToSk < mincentdis) {
						mincentdis = curdisToSk;
						mincent = i;
					}
					dis = std::min(dis, curdis);
					//if (collide && between) {
					if (curdis < l && ratefrom[i]<rate && rate<rateto[i]) {
						voxels[idx].touchid.insert(i);
					}
					if (curdis == 0) {
						voxels[idx].immo = true;
					}
				}
				voxels[idx].hash.assign(bits);
				
				int lev = ((int)dis / (int)far); if (lev >= levlim)lev = levlim - 1;
				if (mincent == -1)lev = levlim - 1;
				if (lev > 0) voxels[idx].tag = Tag(-1, true, lev);
				else voxels[idx].tag = Tag(mincent, bits[mincent], lev);
				
				voxels[idx].setXYZ(ix, iy, iz);
				voxels[idx].dis = dis;
			}
		}
	}
	if (printdebug) { printf("cal vox  - "); toc(); }
	//*******************************************************************
	if(manhmode){
		tic();
		#pragma omp parallel for
		for (int ix = 0; ix < nx; ix++) {
			for (int iy = 0; iy < ny; iy++) {
				for (int iz = 0; iz < nz; iz++) {
					Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
					int idx = ix * ny * nz + iy * nz + iz;
					int itx = ix / stx, ity = iy / sty, itz = iz / stz;
					int idxt = itx * nty * ntz + ity * ntz + itz;
					std::vector<bool> bits = std::vector<bool>(topo.knifes.size());
					for (auto edge : topo.edges) {
						Vector3 v1 = topo.vertices[edge.ia];
						Vector3 v2 = topo.vertices[edge.ib];
						Vector3 vmax = Vector3(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
						Vector3 vmin = Vector3(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
						voxels[idx].hash4.addHash(pos.x - vmax.x < manh);
						voxels[idx].hash4.addHash(pos.x - vmin.x > -manh);
						voxels[idx].hash4.addHash(pos.y - vmax.y < manh);
						voxels[idx].hash4.addHash(pos.y - vmin.y > -manh);
						voxels[idx].hash4.addHash(pos.z - vmax.z < manh);
						voxels[idx].hash4.addHash(pos.z - vmin.z > -manh);
					}
				}
			}
		}
		if (printdebug) { printf("cal manh - "); toc(); }
	}
	//*******************************************************************
	tic();
	//cal grid touch
	std::vector<std::set<int>> tis(ntx * nty * ntz);
	std::vector<std::vector<bool>> tbits2(ntx * nty * ntz, std::vector<bool>(topo.knifes.size(), false));
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				int idx = ix * ny * nz + iy * nz + iz;
				int itx = ix / stx, ity = iy / sty, itz = iz / stz;
				int idxt = itx * nty * ntz + ity * ntz + itz;
				for (auto e : voxels[idx].touchid) {
					tis[idxt].insert(e);
					tbits2[idxt][e] = true;
				}
				voxels[idx].setTXYZ(itx, ity, itz);
				voxels[idx].it = idxt;
			}
		}
	}
	//cal manhatan touch
	std::map<Hash, std::vector<bool>> cubeCollideMap;
	if (manhmode) {
		for (int ix = 0; ix < nx; ix++) for (int iy = 0; iy < ny; iy++) for (int iz = 0; iz < nz; iz++) {
			Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
			int idx = ix * ny * nz + iy * nz + iz;
			Hash hash4 = voxels[idx].hash4;
			if (manhCubeMap.count(hash4) == 0) {
				std::vector<bool> bits4 = hash4.getBits();
				Vector3 _ld = ld;
				Vector3 _ru = ru;
				for (int k = 0; k < topo.edgenum; k++) {
					TopoEdge edge = topo.edges[k];
					Vector3 v1 = topo.vertices[edge.ia];
					Vector3 v2 = topo.vertices[edge.ib];
					Vector3 vmax = Vector3(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z)) + Vector3(manh, manh, manh);
					Vector3 vmin = Vector3(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z)) - Vector3(manh, manh, manh);
					if (bits4[k * 6 + 0])_ru.x = std::min(_ru.x, vmax.x);
					else _ld.x = std::max(_ld.x, vmax.x);
					if (bits4[k * 6 + 1])_ld.x = std::max(_ld.x, vmin.x);
					else _ru.x = std::min(_ru.x, vmin.x);
					if (bits4[k * 6 + 2])_ru.y = std::min(_ru.y, vmax.y);
					else _ld.y = std::max(_ld.y, vmax.y);
					if (bits4[k * 6 + 3])_ld.y = std::max(_ld.y, vmin.y);
					else _ru.y = std::min(_ru.y, vmin.y);
					if (bits4[k * 6 + 4])_ru.z = std::min(_ru.z, vmax.z);
					else _ld.z = std::max(_ld.z, vmax.z);
					if (bits4[k * 6 + 5])_ld.z = std::max(_ld.z, vmin.z);
					else _ru.z = std::min(_ru.z, vmin.z);
				}
				Cube manhcube = Cube(_ld, _ru);
				manhCubeMap[hash4] = manhcube;
				std::vector<bool> cols;
				for (auto cap : topo.caps) {
					cols.push_back(cap.collide(manhcube.cent, manhcube.size.length() / 2));
				}
				cubeCollideMap[hash4] = cols;
			}
		}
	}
	//mod hash
	#pragma omp parallel for
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
				int idx = ix * ny * nz + iy * nz + iz;
				int itx = ix / stx, ity = iy / sty, itz = iz / stz;
				int idxt = itx * nty * ntz + ity * ntz + itz;
				if(!manhmode)voxels[idx].hash2.assign(tbits2[idxt]);
				else voxels[idx].hash2.assign(cubeCollideMap[voxels[idx].hash4]);
				voxels[idx].hash = voxels[idx].hash & voxels[idx].hash2;
				if (tis[idxt].size() == 0)voxels[idx].tag.tlev = 1;
				else voxels[idx].tag.tlev = 0;
			}
		}
	}
	if (printdebug) { printf("mod hash - "); toc(); }
	//*******************************************************************
	supervoxels = std::vector<std::vector<int>>(ntx * nty * ntz);
	supervoxeltouch = std::vector<bool>(ntx * nty * ntz, false);
	for (int idxt = 0; idxt < supervoxeltouch.size(); idxt++) supervoxeltouch[idxt] = (tis[idxt].size() == 0);
	//************************************************************************************
	//cal piece dif in a supervoxel 
	tic();
	std::vector<std::map<Hash, int>> superpieces(ntx * nty * ntz);
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
				int idx = ix * ny * nz + iy * nz + iz;
				int itx = ix / stx, ity = iy / sty, itz = iz / stz;
				int idxt = itx * nty * ntz + ity * ntz + itz;
				Hash hash = voxels[idx].hash & voxels[idx].hash2;
				if (superpieces[idxt].count(hash) == 0)superpieces[idxt][hash]=0;
				superpieces[idxt][hash]++;
			}
		}
	}
	if (printdebug) { printf("cal sup  - "); toc(); }
	//************************************************************************************
	Group group; group.initAva();topo.avalist = group.avalist;
	//************************************************************************************
	int cal = 0;
	for (int i = 0; i < superpieces.size(); i++) {
		int min = INT_MAX, max = INT_MIN;
		for (std::pair<Hash, int> it : superpieces[i]) {
			min = std::min(min, it.second);
			max = std::max(max, it.second);
		}
		cal += max-min;
	}
	return cal;
}

void Utility::genPiece_voxel_bfs() {
	std::set<int> voxels_set;
	for (int i = 0; i < voxelsize; i++) {
		if(voxels[i].dis>0)if(!nearmode || voxels[i].tag.tlev == 0)voxels_set.insert(i);
	}
	while (true) {
		int pn = pieces.size();
		Piece piece;
		int flag = 0;
		std::vector<int> queue;
		int root = -1;
		for (auto i : voxels_set) {
			if (!voxels[i].exist) {
				root = i;
				break;
			}
		}
		if (root == -1)break;
		queue.push_back(root);
		voxels[root].exist = true;
		while (flag < queue.size()) {
			int qn = queue.size();
			for (; flag < qn; flag++) {
				int cur = queue[flag];
				piece.voxelsi.push_back(cur);
				piece.volume++;
				voxels[cur].belong = pn;
				voxels_set.erase(cur);
				int ix, iy, iz;
				voxels[cur].getXYZ(ix, iy, iz);
				int nei[6] = { -1, -1, -1, -1, -1, -1 };
				if (ix - 1 >= 0) nei[0] = ((ix - 1) * ny * nz + (iy) * nz + (iz));
				if (ix + 1 < nx) nei[1] = ((ix + 1) * ny * nz + (iy) * nz + (iz));
				if (iy - 1 >= 0) nei[2] = ((ix) * ny * nz + (iy - 1) * nz + (iz));
				if (iy + 1 < ny) nei[3] = ((ix) * ny * nz + (iy + 1) * nz + (iz));
				if (iz - 1 >= 0) nei[4] = ((ix) * ny * nz + (iy) * nz + (iz - 1));
				if (iz + 1 < nz) nei[5] = ((ix) * ny * nz + (iy) * nz + (iz + 1));
				for (int i = 0; i < 6; i++) {
					if (nei[i] == -1)continue;
					int tar = nei[i];
					if (voxels[tar].exist)continue;
					if (voxels[tar].dis == 0)continue;
					if (nearmode && voxels[tar].tag.tlev != 0)continue;
					if(tagmode)if (voxels[cur].tag == voxels[tar].tag) {
						queue.push_back(tar);
						voxels[tar].exist = true;
					}
					if (!tagmode)if (voxels[cur].hash == voxels[tar].hash && voxels[cur].hash2 == voxels[tar].hash2)
					if(voxels[cur].tequal(voxels[tar])){
						queue.push_back(tar);
						voxels[tar].exist = true;
					}
				}
			}
		}
		//******************************************************
		std::set<TouchInfo> tis;
		for (int i = 0; i < 6; i++)piece.hash.addHash(true);
		for (auto idx : piece.voxelsi) {
			Vector3 pos = voxels[idx].pos;
			for (auto e : voxels[idx].touchid) {
				Vector3 norm = topo.knifes[e].normal.normalize();
				float dis = topo.knifes[e].distanceToPoint(pos);
				Vector3 dir = dis > 0 ? norm : norm * -1;
				if(std::abs(dis)>0.1f)tis.insert(TouchInfo(e, dir));
			}
			piece.hash = piece.hash & voxels[idx].hash3;
		}
		for (TouchInfo ti : tis) {
			piece.touchinfos.push_back(ti);
		}
		piece.id = pn;
		topo.calTouchBound(piece, false);
		if (piece.boundids.size() == topo.avalist.size())
			topo.calTouchBound(piece, true);
		piece.it = voxels[piece.voxelsi[0]].it;
		pieces.push_back(piece);
		supervoxels[piece.it].push_back(pn);
		if(pn % 100==0)std::cout << "Generate Piece" << " " << pn << std::endl;
	}
}

void Utility::genPiece_voxel() {
	//************************************************************************************
	tic(); 
	voxelDirSeen();
	if (printdebug) { printf("cal seen - "); toc(); }
	//genPiece_voxel_bfs();return;
	//************************************************************************************
	bool onlymode = true;
	std::vector<std::map<Hash, int>> hashmaps(ntx*nty*ntz);
	for (int ix = 0; ix < nx; ix++) for (int iy = 0; iy < ny; iy++) for (int iz = 0; iz < nz; iz++) {
		int idx = ix * ny * nz + iy * nz + iz;
		int itx = ix / stx, ity = iy / sty, itz = iz / stz;
		int idxt = itx * nty * ntz + ity * ntz + itz;
		Hash hash = voxels[idx].hash;
		if(manhmode)for (auto b : voxels[idx].hash4.getBits()) hash.addHash(b);
		if (hashmaps[idxt].count(hash) == 0) {
			Piece piece;
			piece.id = pieces.size();
			hashmaps[idxt][hash] = piece.id;
			//pieces.push_back(Piece());
			pieces.push_back(piece);
		}
		int p = hashmaps[idxt][hash];
		voxels[idx].belong = p;
		voxels[idx].exist = true;
		pieces[p].voxelsi.push_back(idx);
		pieces[p].volume++;
	}

	for (int p = 0; p < pieces.size(); p++) {
		std::set<TouchInfo> tis;
		for (int i = 0; i < 6; i++)pieces[p].hash.addHash(true);
		for (auto idx : pieces[p].voxelsi) {
			Vector3 pos = voxels[idx].pos;
			for (auto e : voxels[idx].touchid) {
				Vector3 norm = topo.knifes[e].normal.normalize();
				Vector3 dir = topo.knifes[e].distanceToPoint(pos) > 0 ? norm : norm * -1;
				tis.insert(TouchInfo(e, dir));
			}
			pieces[p].hash = pieces[p].hash & voxels[idx].hash3;
		}
		for (TouchInfo ti : tis) {
			pieces[p].touchinfos.push_back(ti);
		}
		pieces[p].id = p;
		pieces[p].it = voxels[pieces[p].voxelsi[0]].it;
		supervoxels[pieces[p].it].push_back(p);
	}
	if (printdebug) { printf("Generate Piece %d\n", pieces.size()); }
	tic();
	#pragma omp parallel for
	for (int p = 0; p < pieces.size(); p++) {
		topo.calTouchBound(pieces[p], false);
		if (pieces[p].boundids.size() == topo.avalist.size())
			topo.calTouchBound(pieces[p], true);
	}
	if (printdebug) { printf("cal bound- "); toc(); }
	//************************************************************************************
	for(int i=0;i<pieces.size();i++){
		pieces[i].max = Vector3(FLT_MIN, FLT_MIN, FLT_MIN), pieces[i].min = Vector3(FLT_MAX, FLT_MAX, FLT_MAX);
		for (auto v : pieces[i].voxelsi) {
			pieces[i].max.x = std::max(pieces[i].max.x, voxels[v].pos.x);
			pieces[i].max.y = std::max(pieces[i].max.y, voxels[v].pos.y);
			pieces[i].max.z = std::max(pieces[i].max.z, voxels[v].pos.z);
			pieces[i].min.x = std::min(pieces[i].min.x, voxels[v].pos.x);
			pieces[i].min.y = std::min(pieces[i].min.y, voxels[v].pos.y);
			pieces[i].min.z = std::min(pieces[i].min.z, voxels[v].pos.z);
		}
	}
	if (printdebug) { printf("cal max- "); toc(); }
	//************************************************************************************
	/*
	int mapcnt = 0;
	std::map<Tag, int> tagmap;
	for (int i = 0; i < voxelsize; i++) if (tagmap.count(voxels[i].tag) == 0)tagmap[voxels[i].tag] = mapcnt++;
	pieces = std::vector<Piece>(tagmap.size());
	for (int i = 0; i < voxelsize; i++) {
		int h = tagmap[voxels[i].tag];
		pieces[h].voxelsi.push_back(i);
		pieces[h].volume++;
		voxels[i].exist = true;
		voxels[i].belong = h;
	}
	for (int i = 0; i < pieces.size(); i++) {
		int idx = pieces[i].voxelsi[0];
		int e = voxels[idx].tag.e;
		Vector3 norm = topo.knifes[e].normal.normalize();
		Vector3 dir = voxels[idx].tag.side ? norm : norm*-1;
		pieces[i].touchinfos.push_back(TouchInfo(e, dir));
		pieces[i].id = i;
		//std::set<TouchInfo> tis;
		//float maxd = 0;
		//for (int j = 0; j < pieces[i].voxelsi.size(); j++) {
		//	int idx = pieces[i].voxelsi[j];
		//	for (auto t : voxels[idx].touchid) {
		//		Vector3 norm = topo.knifes[t].normal.normalize();
		//		Vector3 dir = voxels[idx].hash.getBit(t) ? norm : norm*-1;
		//		float d = abs(topo.knifes[t].distanceToPoint(voxels[idx].pos));
		//		if (d > maxd) {
		//			tis.insert(TouchInfo(t, dir));//apply last dir
		//			maxd = d;
		//		}
		//	}
		//}
		//for (auto t : tis) {
		//	pieces[i].touchinfos.push_back(t);
		//}
		//pieces[i].id = i;
	}
	*/
	//************************************************************************************
	/* old version
	std::map<Hash, Piece> hashMap;
	for (int i = 0; i < voxelsize; i++) {
		int ix = i / (ny * nz);
		int iy = i / nz % ny;
		int iz = i % nz;
		Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
		//if (voxels[i].dis < topo.radii*3 && voxels[i].dis > topo.radii) {
		if (voxels[i].dis <= far) {
			if (hashMap.count(voxels[i].hash) == 0) {
				hashMap[voxels[i].hash] = Piece(Mesh(), voxels[i].hash);//insert
			}
			hashMap[voxels[i].hash].voxelsi.push_back(i);
			hashMap[voxels[i].hash].volume++;
			voxels[i].exist = true;
		}
	}
	std::map<Hash, Piece>::iterator it;
	for (it = hashMap.begin(); it != hashMap.end(); it++) {
		pieces.push_back(it->second);
	}
	toc();
	//************************************************************************************
	#pragma omp parallel for
	for (int i = 0; i < pieces.size(); i++) {
		for (int j = 0; j < pieces[i].voxelsi.size(); j++) {
			voxels[pieces[i].voxelsi[j]].belong = i;
		}
		pieces[i].id = i;
	}
	toc();
	*/
}

void Utility::genPiece() {
	pieces.push_back(Piece(Mesh().genCube(ld, ru)));
	int n = topo.knifes.size();
	for (int j = 0; j< n; j++) {
		std::cout << "                    " << std::endl;
		std::cout << "Apply Knife " << j + 1 << "/" << n << std::endl;
		std::vector<Piece> pieces_next = std::vector<Piece>();
		#pragma omp parallel for
		for (int i = 0; i < pieces.size(); i++) {
			Piece piece = pieces[i];
			Piece piece_ = Piece(piece.mesh.slice(topo.knifes[j].center, topo.knifes[j].normal), piece.hash);
			if (true) {
				if (piece.mesh.vertices.size() > 0) {
					piece.hash.addHash(true);
					#pragma omp critical
					pieces_next.push_back(piece);
				}
				if (piece_.mesh.vertices.size() > 0) {
					piece_.hash.addHash(false);
					#pragma omp critical
					pieces_next.push_back(piece_);
				}
			}
		}

		pieces.swap(pieces_next);
		std::cout << "Pieces Num : " << pieces.size() << std::endl;
		std::cout << "____________________" << std::endl;
	}
	for (int i = 0; i < pieces.size(); i++) {
		pieces[i].id = i;
	}
}

void Utility::initGroup() {//assume piece.id == group.id
	Group group;group.initAva();
	groups = std::vector<Group>(pieces.size());
	groupIdxMap = std::vector<int>(groups.size());
	idexist = std::vector<bool>(groups.size(), true);
	groupidcnt = pieces.size();
	int toremovenum = 0;
	//#pragma omp parallel for
	for (int i = 0; i < pieces.size(); i++) {
		if (i % 100 == 0)printf("%d / %d\n", i, pieces.size());
		//groups[i] = group; //"no more initava"
		groups[i].avaset = group.avaset;
		//topo.calTouch(pieces[i]);
		appendPiece(groups[i], pieces[i]);
		groups[i].id = pieces[i].id; // == i
		groups[i].initAvaLev(topo.maxcol);
		pieces[i].belong = groups[i].id;
		groupIdxMap[groups[i].id] = i;//groupIdxMap.insert(std::pair<int, int>(groups[i].id, i));
		initGroupBorde(groups[i]);
		if (colmode) {
			std::set<int> ids; ids.insert(i); calFarContact(groups[i], ids);
		}
		//****************************
		int tar = i;
		if (groups[tar].volume < 0) {
			if(printdebug)printf("Small Volume : %f\n", groups[tar].volume);
			for (int i = 0; i < groups[tar].pieces.size(); i++) {
				for (int j = 0; j < pieces[groups[tar].pieces[i]].voxelsi.size(); j++) {
					int idx = pieces[groups[tar].pieces[i]].voxelsi[j];
					voxels[idx].removed = true;
				}
				pieces[groups[tar].pieces[i]].removed = true;
			}
			idexist[groups[tar].id] = false;
			groups[tar].removed = true;
			if (!checkThin(tar)) {
				removeGroup(tar);
				toremovenum++;
			}
		}
		//****************************
	}
	if (printdebug)printf("InitGroup %d, Too Thin : %d \n", groups.size(), toremovenum);
}

void Utility::initGroupBorde(Group &group) {
	group.borde = std::vector<int>();
	for (int i = 0; i < group.pieces.size(); i++) {
		for (int j = 0; j < pieces[group.pieces[i]].voxelsi.size(); j++) {
			int idx = pieces[group.pieces[i]].voxelsi[j];
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			int ix, iy, iz; voxels[idx].getXYZ(ix, iy, iz);
			if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			int samecnt = 0;
			for (int j = 0; j < 6; j++) {//assume piece id is equal to group id
				if (nei[j] == -1) continue;
				if (!voxels[nei[j]].exist)continue;
				if (!voxels[nei[j]].removed)continue;
				if (voxels[idx].belong == voxels[nei[j]].belong)samecnt++;
			}
			if (samecnt < 6)group.borde.push_back(idx);
		}
	}
	group.bordenum = group.borde.size();
}

void Utility::purneGroupBorde(Group &group) {
	int flag = 0;
	for (int i = 0; i < group.bordenum; i++) {
		int idx = group.borde[i];
		int ix, iy, iz; voxels[idx].getXYZ(ix, iy, iz);
		int nei[6] = { -1, -1, -1, -1, -1, -1 };
		if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
		if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
		if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
		if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
		if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
		if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
		int samecnt = 0;
		for (int j = 0; j < 6; j++) {//assume piece id is equal to group id
			if (nei[j] == -1) continue;
			if (!voxels[nei[j]].exist)continue;
			if (!voxels[nei[j]].removed)continue;
			if (voxels[idx].belong == voxels[nei[j]].belong)samecnt++;
		}
		if (samecnt < 6)group.borde[flag++] = idx;
	}
	group.bordenum = flag;
}

void Utility::calTouch(int p) {//assume all knifes and only touch
	std::set<TouchInfo> touchinfos;
	for (int i = 0; i < pieces[p].voxelsi.size();i++) {
		int idx = pieces[p].voxelsi[i];
		Vector3 pos = voxels[idx].pos;
		for (int j = 0; j < topo.edgenum; j++) {
			Vector3 norm = topo.knifes[j].normal.normalize();
			Vector3 dir = voxels[idx].hash.getBit(j)? norm : norm*-1;
			bool between;
			bool collide = topo.caps[j].collide(pos, -l, between);
			if (collide && between) {
				float l1 = (topo.caps[j].proj(pos) - topo.caps[j].p1).length();
				float l = (topo.caps[j].p2 - topo.caps[j].p1).length();
				if ((l1 / l) > 0.3 && (l1 / l) < 0.7) {
					touchinfos.insert(TouchInfo(j, dir));
				}
			}
		}
	}
	//std::cout << touchinfos.size() << std::endl;
	pieces[p].touchinfos = std::vector<TouchInfo>();
	for (auto ti : touchinfos) {
		pieces[p].touchinfos.push_back(ti);
	}
}

void Utility::initLink() {
	for (int i = 0; i < groups.size(); i++) {
		std::vector<int> forces = std::vector<int>(groups.size());
		#pragma omp parallel for
		for (int j = i + 1; j < groups.size(); j++) {
			forces[j] = 0;
			for (int n = 0; n < groups[i].pieces.size(); n++) {
				for (int m = 0; m < groups[j].pieces.size(); m++) {
					if (pieces[groups[i].pieces[n]].isNei(pieces[groups[j].pieces[m]])) {
						forces[j]++;
					}
				}
			}
		}
		for (int j = i + 1; j < groups.size(); j++) {
			if (forces[j] > 0) {
				groups[i].neighbor.insert(groups[j].id);
				groups[j].neighbor.insert(groups[i].id);
				groupLink.push(GroupLink(groups[i].id, groups[j].id, calWorth(groups[i], groups[j])));
			}
		}
	}
}

void Utility::initLink_voxel() {//assume group.id = piece.id
	tic();
	//#pragma omp parallel for
	for (int i = 0; i < voxelsize; i++) {
		int ix = i / (ny * nz);
		int iy = i / nz % ny;
		int iz = i % nz;
		int nei[6] = { -1, -1, -1, -1, -1, -1 };
		if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy) * nz + (iz));
		if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy) * nz + (iz));
		if (iy - 1 >= 0)nei[2] = ((ix) * ny * nz + (iy - 1) * nz + (iz));
		if (iy + 1 < ny)nei[3] = ((ix) * ny * nz + (iy + 1) * nz + (iz));
		if (iz - 1 >= 0)nei[4] = ((ix) * ny * nz + (iy) * nz + (iz - 1));
		if (iz + 1 < nz)nei[5] = ((ix) * ny * nz + (iy) * nz + (iz + 1));
		int samecnt = 0;
		for (int j = 0; j < 6; j++) {//assume piece id is equal to group id
			if (nei[j] != -1) {
				if (!voxels[i].exist || !voxels[nei[j]].exist)continue;
				int a = voxels[i].belong;
				int b = voxels[nei[j]].belong;
				if (a != b) {
					if (pieces[a].contactmap.count(b) == 0)pieces[a].contactmap[b] = std::vector<int>(6, 0);
					pieces[a].contactmap[b][j]++;
					if (j == 4) pieces[a].downvoxels.push_back(i);
					if (j == 5) pieces[a].upvoxels.push_back(i);
				}
			} else {
				int a = voxels[i].belong;
				if (j == 4) pieces[a].downvoxels.push_back(i);
				if (j == 5) pieces[a].upvoxels.push_back(i);
			}
		}
	}
	if (printdebug) { std::cout << "cal pcnt - "; toc(); }
	//***********************************************************************
	if (checkcavmode) {
		#pragma omp parallel for
		for (int i = 0; i < groups.size(); i++) {
			for (auto p : groups[i].pieces) {
				checkCav(groups[i], pieces[p]);
			}
			collectCav(groups[i]);
		}
	}
	//***********************************************************************
	tic();
	//#pragma omp parallel for
	for(int i=0;i<pieces.size();i++){
		for (auto jj : pieces[i].contactmap) {
			int j = jj.first;
			std::vector<int> contact= jj.second;
			/*
			int contactface = 0;
			for (int k = 0; k < 6; k++) {
				contactface += piececontact[i][j][k];
			}
			if (contactface > 0) {
				groups[i].contactmap[j] = piececontact[i][j];
			}
			*/
			//for (int k = 0; k < 6; k++)printf("%d ", contact[k]); printf("\n");
			pieces[i].neighbor.insert(pieces[j].id);
			groups[i].contactmap = pieces[i].contactmap;
			groups[i].neighbor.insert(groups[j].id);
		}
	}
	if (printdebug) { std::cout << "cal neig - "; toc(); }
	//***********************************************************************
	tic();
	for (int i = 0; i < groups.size(); i++) {
		for (auto jj : groups[i].neighbor) {
			int j = groupIdxMap[jj];
			if (i >= j)continue;
			groupLink.push(GroupLink(groups[i].id, groups[j].id, calWorth(groups[i], groups[j])));
		}
	}
	if (printdebug) { std::cout << "cal gplk - "; toc(); }
}

void Utility::MergeGroup(Group& group1, Group& group2) {
	Group group = group1;
	for (int i = 0; i < group2.pieces.size(); i++) {
		appendPiece(group, pieces[group2.pieces[i]]);
	}
	//checkCav(group, group2);

	//*********************************************//voxels & piece belong
	for (int i = 0; i < group1.pieces.size(); i++)
		for (int j = 0; j < pieces[group1.pieces[i]].voxelsi.size(); j++)
			voxels[pieces[group1.pieces[i]].voxelsi[j]].belong = groupidcnt;
	for (int i = 0; i < group2.pieces.size(); i++)
		for (int j = 0; j < pieces[group2.pieces[i]].voxelsi.size(); j++)
			voxels[pieces[group2.pieces[i]].voxelsi[j]].belong = groupidcnt;
	group.id = groupidcnt++;
	for (auto p : group.pieces)pieces[p].belong = group.id;
	//*********************************************//group borde
	/*
	group.borde = std::vector<int>(group1.bordenum + group2.bordenum);
	for (int i = 0; i < group1.bordenum; i++)group.borde[i] = group1.borde[i];
	for (int i = 0; i < group2.bordenum; i++)group.borde[i+group1.bordenum] = group2.borde[i];
	group.bordenum = group.borde.size();
	purneGroupBorde(group);
	*/
	//*********************************************//group neighbor
	idexist.push_back(true);
	idexist[group1.id] = idexist[group2.id] = false;
	group.neighbor.swap(std::set<int>());
	for (auto id : group1.neighbor)if (id != group1.id && id != group2.id && idexist[id]) { group.neighbor.insert(id);}
	for (auto id : group2.neighbor)if (id != group1.id && id != group2.id && idexist[id]) { group.neighbor.insert(id);}
	//*********************************************//group contact map
	group.contactmap.clear(); //"don't copy in operator=" is useless
	Group* grouplist[2] = { &group1, &group2 };
	for (int n = 0; n < 2; n++) {
		for (auto iter : grouplist[n]->contactmap) {
			if (!idexist[iter.first])continue;//if (iter.first == group1.id || iter.first == group2.id)continue;
			if (group.contactmap.count(iter.first) != 0) {
				std::vector<int> contact;
				contact.swap(group.contactmap[iter.first]);
				for (int i = 0; i < 6; i++) {
					contact[i] += iter.second[i];
				}
				group.contactmap[iter.first].swap(contact);
			}
			else {
				group.contactmap[iter.first] = iter.second;
			}
		}
	}
	for (auto id : group.neighbor) {
		std::vector<int> contact(6, 0);
		for (int n = 0; n < 2; n++) {
			if (groups[id].contactmap.count(grouplist[n]->id) > 0) {
				std::vector<int> *temp = &groups[id].contactmap[grouplist[n]->id];
				for (int i = 0; i < 6; i++) {
					contact[i] += (*temp)[i];
				}
			}
		}
		groups[id].contactmap[group.id].swap(contact);
	}
	
	//*******************************//use idexist to check, so no need
	//std::vector<Group> groups_new;
	//for (int i = 0; i < groups.size(); i++) {
	//	if (groups[i].id != group1.id && groups[i].id != group2.id) {
	//		//groups[i].neighbor.erase(group1.id);//用idexist加速
	//		//groups[i].neighbor.erase(group2.id);
	//		groups_new.push_back(groups[i]);
	//		groupIdxMap[groups[i].id] = (int)(groups_new.size() - 1);
	//	}
	//}
	//groups.swap(groups_new);
	//groupIdxMap.erase(group1.id);//use idexist to check
	//groupIdxMap.erase(group2.id);//now group1,2 is not avaliable now

	//*********************************************//group coldis
	if(colmode)appendGroup(group, group2);
	//*********************************************//others
	group.trigger = true;
	groups.push_back(group);
	groupIdxMap.push_back((int)(groups.size() - 1));//groupIdxMap.insert(std::pair<int, int>(group.id, (int)(groups.size() - 1)));
	//purneVoxel(group);
	
	if(printdebug_l2)std::cout << "New Neighbor number : " << group.neighbor.size() << std::endl;
	for (auto id : group.neighbor) {
		groups[groupIdxMap[id]].neighbor.insert(group.id);
		Worth worth = calWorth(group, groups[groupIdxMap[id]]);
		float volume = (group.volume + groups[groupIdxMap[id]].volume);
		groupLink.push(GroupLink(group.id, id, worth));
	}
}

void Utility::checkCav(Group& group, Group& group2) {
	if (group.downplane.size() == 0) group.downplane = std::vector<std::set<int>>(nx * ny, std::set<int>());
	if (group.upplane.size() == 0) group.upplane = std::vector<std::set<int>>(nx * ny, std::set<int>());
	std::set<int> newixy;
	for (auto i : group2.downvoxels) {
		int ix, iy, iz; voxels[i].getXYZ(ix, iy, iz);
		int ixy = ix * ny + iy;
		newixy.insert(ixy);
		group.downplane[ixy].insert(iz);
	}
	for (auto i : group2.upvoxels) {
		int ix, iy, iz; voxels[i].getXYZ(ix, iy, iz);
		int ixy = ix * ny + iy;
		newixy.insert(ixy);
		group.upplane[ixy].insert(iz);
	}
	for (auto ixy : newixy) {
		std::set<int> downremove;
		std::set<int> upremove;
		for (auto a : group.downplane[ixy]) {
			for (auto b : group.upplane[ixy]) {
				if (a == b + 1) {
					downremove.insert(a);
					upremove.insert(b);
				}
			}
		}
		for (auto r : downremove) group.downplane[ixy].erase(r);
		for (auto r : upremove) group.upplane[ixy].erase(r);
	}
	group.downvoxels.clear();
	group.upvoxels.clear();
	int upcnt = 0;
	int downcnt = 0;
	for (auto ixy : newixy) {
		if (group.downplane[ixy].size() > 1)downcnt++;
		if (group.upplane[ixy].size() > 1)upcnt++;
	}
	if(downcnt > colthre || upcnt > colthre)group.zcav = true;
	collectCav(group);
}

void Utility::collectCav(Group& group) {
	if (!quickcheckcav) {
		for (int ix = 0; ix < nx; ix++)for (int iy = 0; iy < ny; iy++) {
			int ixy = ix * ny + iy;
			for (auto iz : group.downplane[ixy]) group.downvoxels.push_back(ixy * nz + iz);
			for (auto iz : group.upplane[ixy]) group.upvoxels.push_back(ixy * nz + iz);
		}
	}
	else {
	
	}
}

void Utility::checkCav(Group& group) {
	if (!quickcheckcav) {

	} else{
		for (int itxy = 0; itxy < ntx*nty; itxy++) {
			int upnum = 0;
			int downnum = 0;
			for (auto itzmap : group.ztmap[itxy]) {
				int it = itxy * ntz + itzmap.first;
				if (itzmap.first + 1 >= ntz)upnum++;
				else if (itzmap.first - 1 < 0)downnum++;
				else {
					int lim = 50;
					for (auto p : group.supervoxels[it]) {
						int it_ = itxy * ntz + itzmap.first + 1;
						bool isnei = false;
						if (group.supervoxels.count(it_) > 0) {
							for (auto p_ : group.supervoxels[it_]) {
								if (pieces[p].neighbor.count(p_) > 0) {
									int cn = 0;
									for (auto n : pieces[p].contactmap[p_]) {
										cn += n;
									}
									if(cn>lim)isnei = true;
								}
							}
						}
						if (!isnei)upnum++;
						it_ = itxy * ntz + itzmap.first - 1;
						isnei = false;
						if (group.supervoxels.count(it_) > 0) {
							for (auto p_ : group.supervoxels[it_]) {
								if (pieces[p].neighbor.count(p_) > 0) {
									int cn = 0;
									for (auto n : pieces[p].contactmap[p_]) {
										cn += n;
									}
									if (cn>lim)isnei = true;
								}
							}
						}
						if (!isnei)downnum++;
					}
				}
			}				
			if (downnum > 1 || upnum > 1) {
				group.zcav = true;
				return;
			}
		}
		return;//force try above
		//*******************************
		for(auto set : group.ztmap){
			if (set.size() <= 1)continue;
			std::vector<int> exist(ntz, 0);
			for (auto itzmap : set) {
				exist[itzmap.first] = itzmap.second;
			}
			int es = exist.size();
			int upnum = 0;
			int downnum = 0;
			int vlim = 100;
			int cubv = stx*sty*stz;

			int cur = 0;
			bool dec = false;
			for (int i = 0; i < es - 1; i++) {
				if (!dec) {
					if (cur > exist[i])dec = true;
					cur = exist[i];
				}
				else {
					if (cur < exist[i]) {
						dec = false;
						downnum++;
					}
					cur = exist[i];
				}
			}
			cur = 0;
			dec = false;
			for (int i = es - 1; i > 0; i--) {
				if (!dec) {
					if (cur > exist[i])dec = true;
					cur = exist[i];
				}
				else {
					if (cur < exist[i]) {
						dec = false;
						upnum++;
					}
					cur = exist[i];
				}
			}
			if (downnum >= 1 || upnum >= 1) {
				group.zcav = true;
				return;
			}
		}
	}
}

void Utility::checkCav(Group& group, Piece& piece) {
	if (!quickcheckcav) {
		if (group.downplane.size() == 0) group.downplane = std::vector<std::set<int>>(nx * ny, std::set<int>());
		if (group.upplane.size() == 0) group.upplane = std::vector<std::set<int>>(nx * ny, std::set<int>());
		std::set<int> newixy;
		for (auto i : piece.downvoxels) {
			int ix, iy, iz; voxels[i].getXYZ(ix, iy, iz);
			int ixy = ix * ny + iy;
			newixy.insert(ixy);
			group.downplane[ixy].insert(iz);
		}
		for (auto i : piece.upvoxels) {
			int ix, iy, iz; voxels[i].getXYZ(ix, iy, iz);
			int ixy = ix * ny + iy;
			newixy.insert(ixy);
			group.upplane[ixy].insert(iz);
		}
		for (auto ixy : newixy) {
			std::set<int> downremove;
			std::set<int> upremove;
			for (auto a : group.downplane[ixy]) {
				for (auto b : group.upplane[ixy]) {
					if (a == b + 1) {
						downremove.insert(a);
						upremove.insert(b);
					}
				}
			}
			for (auto r : downremove) group.downplane[ixy].erase(r);
			for (auto r : upremove) group.upplane[ixy].erase(r);
		}
		for (auto ixy : newixy) {//Big Problem!, group will be ok while piece may not
			if (group.downplane[ixy].size() > 1 || group.upplane[ixy].size() > 1)group.zcav = true;
			/*
			for (auto a : group.downplane[ixy]) {
			for (auto b : group.upplane[ixy]) {
			if (a > b) {
			group.zcav = true;
			}
			}
			}
			*/
		}
	}
	else {
		if (group.ztmap.size() == 0) group.ztmap = std::vector<std::map<int, int>>(ntx * nty, std::map<int, int>());
		int itx, ity, itz; voxels[piece.voxelsi[0]].getTXYZ(itx, ity, itz);
		int itxy = itx * nty + ity;
		if (group.ztmap[itxy].count(itz) == 0)group.ztmap[itxy][itz] = 0;
		group.ztmap[itxy][itz]+=piece.volume;
		//***************************************
		int it = voxels[piece.voxelsi[0]].it;
		if (group.supervoxels.count(it) == 0)group.supervoxels[it] = std::set<int>();
		group.supervoxels[it].insert(piece.id);
	}
	
}

void Utility::appendPiece(Group& group, Piece& piece){
	group.pieces.push_back(piece.id); 
	group.volume += piece.volume;
	topo.purneAva(group, piece);
	//******************************************************
	group.max.x = std::max(piece.max.x, group.max.x);
	group.max.y = std::max(piece.max.y, group.max.y);
	group.max.z = std::max(piece.max.z, group.max.z);
	group.min.x = std::min(piece.min.x, group.min.x);
	group.min.y = std::min(piece.min.y, group.min.y);
	group.min.z = std::min(piece.min.z, group.min.z);
	//******************************************************
	group.idxts.insert(voxels[piece.voxelsi[0]].it);
	//******************************************************
	if(checkcavmode)checkCav(group, piece);
}

void Utility::appendGroup(Group& group, Group& group2) {
	for (int mode = 0; mode < 6; mode++) {
		for (auto pair : group2.idxtcoldis[mode]) {
			std::vector<int> coldis;
			if (group.idxtcoldis[mode].count(pair.first) == 0) {
				group.idxtcoldis[mode][pair.first] = pair.second;
			}
			else {
				/*
				coldis = group.idxtcoldis[mode][pair.first];
				for (int i = 0; i < coldis.size(); i++) {
					coldis[i] = std::min(coldis[i], pair.second[i]);
				}
				group.idxtcoldis[mode][pair.first].swap(coldis);
				*/
				group.idxtcoldis[mode].erase(pair.first);
			}
		}
	}
}

void Utility::calBound() {
	ld = Vector3(FLT_MAX, FLT_MAX, FLT_MAX);
	ru = Vector3(FLT_MIN, FLT_MIN, FLT_MIN);
	for (int i = 0; i < topo.vertices.size(); i++) {
		ld.x = ld.x < topo.vertices[i].x ? ld.x : topo.vertices[i].x;
		ld.y = ld.y < topo.vertices[i].y ? ld.y : topo.vertices[i].y;
		ld.z = ld.z < topo.vertices[i].z ? ld.z : topo.vertices[i].z;
		ru.x = ru.x > topo.vertices[i].x ? ru.x : topo.vertices[i].x;
		ru.y = ru.y > topo.vertices[i].y ? ru.y : topo.vertices[i].y;
		ru.z = ru.z > topo.vertices[i].z ? ru.z : topo.vertices[i].z;
	}
	float radii = topo.radii;
	if (FILE *file = fopen("boundinfo.txt", "r")) {
		if(printdebug)std::cout << "Read Bound From File......" << std::endl;
		fscanf(file, "%f %f %f", &ld.x, &ld.y, &ld.z);
		fscanf(file, "%f %f %f", &ru.x, &ru.y, &ru.z);
		fscanf(file, "%d %d %d", &stx, &sty, &stz);
		nx = (ru.x - ld.x + 0.1) / l;
		ny = (ru.y - ld.y + 0.1) / l;
		nz = (ru.z - ld.z + 0.1) / l;
		ntx = (nx + stx - 1) / stx;
		nty = (ny + sty - 1) / sty;
		ntz = (nz + stz - 1) / stz;
		fclose(file);
	}
	else {
		ld -= Vector3(radii * 2, radii * 2, radii * 2);
		ru += Vector3(radii * 2, radii * 2, radii * 2);
		nx = (ru.x - ld.x) / l + 1;
		ny = (ru.y - ld.y) / l + 1;
		nz = (ru.z - ld.z) / l + 1;
		ntx = (nx + stx - 1) / stx;
		nty = (ny + sty - 1) / sty;
		ntz = (nz + stz - 1) / stz;
		nx = ntx * stx;
		ny = nty * sty;
		nz = ntz * stz;
		Vector3 tune = Vector3(ld.x + nx*l, ld.y + ny*l, ld.z + nz*l) - ru;
		ru += move;
		ld -= move;
		ru += tune / 2;
		ld -= tune / 2;
	}
	//********************************
	//voxels = std::vector<Voxel>(nx*ny*nz);
	voxelsize = nx*ny*nz;
	voxels = new Voxel[voxelsize];
	tnear = std::vector<std::vector<int>>(ntx * nty * ntz);
	printf("%d voxels", nx*ny*nz);
	//********************************
	topo.maxcol = std::max(nx, std::max(ny, nz));
	if (printdebug) {
		std::cout << ld << std::endl;
		std::cout << ru << std::endl;
	}
	//********************************
	if (manhmode || nogrid) {
		ntx = 1;
		nty = 1;
		ntz = 1;
		stx = nx;
		sty = ny;
		stz = nz;
	}	
	//********************************
	std::vector<bool> isexist(topo.edgenum, false);
	for (int i = 0; i < topo.edgenum; i++) {
		Vector3 v1 = topo.vertices[topo.edges[i].ia];
		Vector3 v2 = topo.vertices[topo.edges[i].ib];
		if ((v1.z > ld.z && v1.z<ru.z) || (v2.z>ld.z && v2.z < ru.z)) {
			isexist[i] = true;
		}
	}
	std::vector<TopoEdge> edges;
	std::vector<Vector3> splitNorm;
	for (int i = 0; i < topo.edgenum; i++) {
		if (isexist[i]) {
			edges.push_back(topo.edges[i]);
			splitNorm.push_back(topo.splitNorm[i]);
		}
	}
	topo.edges = edges;
	topo.splitNorm = splitNorm;
	topo.splitNorm_ori = splitNorm;
	topo.edgenum = edges.size();
	topo.prepareData();
	topo.fixAngleDistance();
	topo.genAllKnife();
	topo.genCapsule();
	//********************************
}

void Utility::optimize() {
	int cnt = 0;
	for (int i = 0; i < groups.size(); i++)if (idexist[groups[i].id])cnt++;
	while (true) {
		if (groupLink.size() == 0)break;
		GroupLink gl = groupLink.top();
		groupLink.pop();
		//std::cout << "Top link energy = " << gl.worth << std::endl;
		if (!idexist[gl.ida]|| !idexist[gl.idb])continue; //if (groupIdxMap.count(gl.ida) == 0 || groupIdxMap.count(gl.idb) == 0)continue;
		if (gl.worth.val[0] == 0)break;
		//std::cout << calWorth(groups[groupIdxMap[gl.ida]], groups[groupIdxMap[gl.idb]]) << std::endl;
		MergeGroup(groups[groupIdxMap[gl.ida]], groups[groupIdxMap[gl.idb]]);
		
		//std::cout << calWorth(groups[groups.size() - 1]) << std::endl;

		cnt--;
		if(printdebug_l2)std::cout << "Group number left : " << cnt << std::endl;
		if (cnt <= 1)break;
	}
}

bool Utility::checkThin(int tar) {
	int minx = INT_MAX;
	int miny = INT_MAX;
	int minz = INT_MAX;
	int maxx = INT_MIN;
	int maxy = INT_MIN;
	int maxz = INT_MIN;
	for (auto p : groups[tar].pieces) {
		for (auto v : pieces[p].voxelsi) {
			minx = std::min(voxels[v].ix, minx);
			miny = std::min(voxels[v].iy, miny);
			minz = std::min(voxels[v].iz, minz);
			maxx = std::max(voxels[v].ix, maxx);
			maxy = std::max(voxels[v].iy, maxy);
			maxz = std::max(voxels[v].iz, maxz);
		}
	}
	int mindif = std::min(maxx - minx, std::min(maxy - miny, maxz - minz));
	if (mindif < 3) { printf("Detect thin : %d - %d\n", tar, mindif); return false; }
	return true;
}

void Utility::removeGroup(int tar) {
	//if(printdebug)std::cout << "---------- Remove " << tar << std::endl;
	//************
	for (int i = 0; i < groups[tar].pieces.size(); i++) {
		for (int j = 0; j < pieces[groups[tar].pieces[i]].voxelsi.size(); j++) {
			int idx = pieces[groups[tar].pieces[i]].voxelsi[j];
			voxels[idx].removed = true;
		}
		pieces[groups[tar].pieces[i]].removed = true;
	}
	idexist[groups[tar].id] = false;
	groups[tar].removed = true;
	//************
	for (int i = 0; i < groups.size(); i++) {
		if(groups[i].idxtcoldis.size()>0)
			groups[i].idxtcoldis.swap(std::vector<std::map<int, std::vector<int>>>(6));
	}
}

void Utility::iterate() {
	while (true) {
		Worth maxmyworth;
		int maxi = -1;
		for (int i = 0; i < groups.size();i++) {
			if (!idexist[groups[i].id])continue;
			if (groups[i].removed)continue;
			Worth worth = calWorth(groups[i]);
			if(worth.val[0] == 0)continue;
			std::set<int> ids; ids.insert(groups[i].id);
			std::vector<int> contact = calContact(groups[i], ids);
			int contactface = 0; for (auto c : contact)if (c > 0)contactface++;
			int contactarea = 0; for (auto c : contact)contactarea += c;
			int minz = INT_MAX;
			int maxz = INT_MIN;
			for (auto p : groups[i].pieces)minz = std::min(voxels[pieces[p].voxelsi[0]].itz, minz);
			for (auto p : groups[i].pieces)maxz = std::max(voxels[pieces[p].voxelsi[0]].itz, maxz);
			//if (contact[5] > 0) contactarea = 0;
			Worth myworth;
			myworth.val.push_back(groups[i].maxDirDis);
			myworth.val.push_back(contactarea);
			//myworth.val.push_back(maxz);
			//myworth.val.push_back(-minz);
			if (maxi == -1 || myworth > maxmyworth) {
				//std::cout << "------------------  " << minz << std::endl;
				maxmyworth = myworth;
				maxi = i;
			}
		}
		int tar = maxi;
		if (tar == -1) break;
		//*******************
		if (colmode) {
			std::set<int> gid; gid.insert(groups[tar].id);
			std::vector<int> coldis = calFarContact(groups[tar], gid);
			topo.applyColdis(groups[tar], coldis);
		}
		//*******************
		if (printdebug)std::cout << "---------- Order Remove " << tar << ", volume : " << maxmyworth.val[1] << std::endl;
		removeGroup(tar);
		if (groups[tar].volume > volumelim && checkThin(tar)) {
			groups_final.push_back(tar);
			finalgroupset.insert(tar);
			purneAvaByNei(tar);
			//printf("%d %d %d %d %d %d\n", groups[tar].avalevel[0], groups[tar].avalevel[1], groups[tar].avalevel[2], groups[tar].avalevel[3], groups[tar].avalevel[4], groups[tar].avalevel[5]);
		}
		//*****************************************************************************
		std::priority_queue<GroupLink> groupLink_;
		int size = groupLink.size();
		while (size>0) {
			GroupLink gl = groupLink.top();
			groupLink.pop(); size--;
			if (idexist[gl.ida] && idexist[gl.idb]) {
				gl.worth = calWorth(groups[gl.ida], groups[gl.idb]);
				groupLink_.push(gl);
			}
		}
		groupLink.swap(groupLink_);
		if(groupLink.size()>0)optimize();
		else break;
	}
}

void Utility::recalAssem() {
	//recover removed
	for (auto g : groups_final) {
		idexist[g] = true;
		groups[g].removed = false;
		for (int i = 0; i < groups[g].pieces.size(); i++) {
			for (int j = 0; j < pieces[groups[g].pieces[i]].voxelsi.size(); j++) {
				int idx = pieces[groups[g].pieces[i]].voxelsi[j];
				voxels[idx].removed = false;
			}
			pieces[groups[g].pieces[i]].removed = false;
		}
		/***/
		groups[g].contactmap = std::map<int, std::vector<int>>();
		groups[g].neighbor = std::set<int>();
		for(auto i : groups[g].pieces){
			for (auto jj : pieces[i].contactmap) {
				int j = jj.first;
				int k = pieces[j].belong;
				if (k == g)continue;
				if (groups[g].contactmap.count(k) == 0) groups[g].contactmap[k] = jj.second;
				else {
					for (int n = 0; n < 6; n++) {
						groups[g].contactmap[k][n]+=jj.second[n];
					}
				}
				groups[g].neighbor.insert(k);
			}
		}
	}
	//*************************************************************
	while (true) {
		Worth worth;
		int mini = -1;

		for (int i = 0; i < groups.size(); i++) {
			if (!idexist[groups[i].id])continue;
			if (groups[i].removed)continue;
			//if (calContact(groups[i])[5] > colthre) continue;
			int curminz = INT_MAX, curmaxz = INT_MIN;
			for (auto p : groups[i].pieces) {
				for (auto v : pieces[p].voxelsi) {
					curminz = std::min(voxels[v].iz, curminz);
					curmaxz = std::max(voxels[v].iz, curmaxz);
				}
			}
			Worth curworth;
			curworth.val.push_back(groups[i].volume > 2000 ? 1 : 0);//???
			curworth.val.push_back(-calContact(groups[i])[5]);
			curworth.val.push_back(curmaxz);
			curworth.val.push_back(curminz);
			curworth.val.push_back(-groups[i].volume);
			if (mini==-1 || worth < curworth) {
				worth = curworth;
				mini = i;
			}
		}
		if (mini == -1) break;
		if (printdebug)std::cout << "---------- Assem " << mini << " zup : " << worth.val[0] << " zdo : " << worth.val[1] << " vol : " << groups[mini].volume << std::endl;
		for (int i = 0; i < groups[mini].pieces.size(); i++) {
			for (int j = 0; j < pieces[groups[mini].pieces[i]].voxelsi.size(); j++) {
				int idx = pieces[groups[mini].pieces[i]].voxelsi[j];
				voxels[idx].removed = true;
			}
			pieces[groups[mini].pieces[i]].removed = true;
		}
		idexist[groups[mini].id] = false;
		groups[mini].removed = true;
		groups_final_assem.push_back(mini);
	}
}

void Utility::noopt() {
	while (true) {
		Worth maxmyworth;
		int maxi = -1;
		for (int i = 0; i < groups.size(); i++) {
			if (!idexist[groups[i].id])continue;
			if (groups[i].removed)continue;
			Worth worth = calWorth(groups[i]);
			if (worth.val[0] == 0)continue;
			Worth myworth;
			myworth.val.push_back(groups[i].maxDirDis);
			myworth.val.push_back(groups[i].volume);
			if (maxi == -1 || myworth > maxmyworth) {
				maxmyworth = myworth;
				maxi = i;
			}
		}
		if (maxi == -1) break;
		if (printdebug)std::cout << "---------- Remove " << maxi << ", volume : " << maxmyworth.val[1] << std::endl;
		//*******************
		if (colmode) {
			std::set<int> gid; gid.insert(groups[maxi].id);
			std::vector<int> coldis = calFarContact(groups[maxi], gid);
			topo.applyColdis(groups[maxi], coldis);
		}
		//*******************
		removeGroup(maxi);
		//*******************
		if (groups[maxi].volume > volumelim && checkThin(maxi)) {
			groups_final.push_back(maxi);
			finalgroupset.insert(maxi);
			purneAvaByNei(maxi);
		}
		if (groups[maxi].volume < 50) {
			//if (printdebug)printf("small : %d\n", groups[maxi].volume);
		}

		//std::cout << "           Ava    " << groups[maxi].avanum << std::endl;std::endl;
	}
	//*********************************************************************
	
	if (nofinalmerge) {
		int nonemovablenum = 0;
		for (int k = 0; k < groups.size(); k++) {
			if (!idexist[groups[k].id])continue;
			if (groups[k].removed)continue;
			if (groups[k].volume < volumelim || !checkThin(k)) {
				removeGroup(k);
				continue;
			}
			if (printdebug)std::cout << "---------- Final Remove " << k << ", volume : " << groups[k].volume << std::endl;
			groups_final.push_back(k);
			finalgroupset.insert(k);
			removeGroup(k);
			nonemovablenum++;
		}
		printf("NonRemovableNum : %d", nonemovablenum);
		return;
	}
	Group group;
	std::vector<int> grouplistidx;
	group.initAva();
	for (int k = 0; k < groups.size(); k++) {
		if (!idexist[groups[k].id])continue;
		if (groups[k].removed)continue;
		idexist[groups[k].id] = false;
		for (int i = 0; i < groups[k].pieces.size(); i++) {
			appendPiece(group, pieces[groups[k].pieces[i]]);
		}
		grouplistidx.push_back(k);
	}
	if (group.pieces.size() > 0) {
		group.id = groupidcnt++;
		idexist.push_back(true);
		groupIdxMap.push_back((int)(groups.size() - 1));
		for (int i = 0; i < group.pieces.size(); i++)
			for (int j = 0; j < pieces[group.pieces[i]].voxelsi.size(); j++)
				voxels[pieces[group.pieces[i]].voxelsi[j]].belong = group.id;
		//*********************************************************************
		//not always on use, not debug yet
		group.contactmap.clear();
		Group** grouplist;
		grouplist = (Group**)malloc(sizeof(Group*) * grouplistidx.size());
		for (int i = 0; i < grouplistidx.size(); i++)grouplist[i] = &groups[grouplistidx[i]];
		for (int n = 0; n < grouplistidx.size(); n++) {
			for (auto iter : grouplist[n]->contactmap) {
				if (!idexist[iter.first] && finalgroupset.count(iter.first)==0)continue;//if (iter.first == group1.id || iter.first == group2.id)continue;
				group.neighbor.insert(iter.first);
				if (group.contactmap.count(iter.first) != 0) {
					std::vector<int> contact;
					contact.swap(group.contactmap[iter.first]);
					for (int i = 0; i < 6; i++) {
						contact[i] += iter.second[i];
					}
					group.contactmap[iter.first].swap(contact);
				}
				else {
					group.contactmap[iter.first] = iter.second;
				}
			}
		}
		//*********************************************************************
		for (auto id : group.neighbor) {
			std::vector<int> contact(6, 0);
			for (int n = 0; n < grouplistidx.size(); n++) {
				if (groups[id].contactmap.count(grouplist[n]->id) > 0) {
					std::vector<int> *temp = &groups[id].contactmap[grouplist[n]->id];
					for (int i = 0; i < 6; i++) {
						contact[i] += (*temp)[i];
					}
				}
			}
			groups[id].contactmap[group.id].swap(contact);
		}
		//*********************************************************************
		groups.push_back(group);
		groups_final.push_back(groups.size() - 1);
		finalgroupset.insert(groups.size() - 1);
	}
}

void Utility::purneAvaByNei(int i) {
	std::set<int> ids; ids.insert(groups[i].id);
	topo.boundAva(groups[i], genBound(calContact(groups[i], ids)));
	//topo.boundFarAva(groups[i], genBound(calFarContact(groups[i], ids)));
	//topo.boundAva(groups[i], genBound(groups[i], ids));
}

std::vector<Vector3> Utility::genBound(Group & group, std::set<int> ids){
	std::vector<Vector3> bound;
	for (int k = 0; k < group.pieces.size(); k++) {
		for (int j = 0; j < pieces[group.pieces[k]].voxelsi.size(); j++) {
			int idx = pieces[group.pieces[k]].voxelsi[j];
			int ix, iy, iz;
			voxels[idx].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1) * nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1) * nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			for (int i = 0; i < 6; i++) {
				if (nei[i] == -1)continue;
				if (!voxels[nei[i]].exist)continue;
				if (voxels[nei[i]].removed)continue;
				if (ids.count(voxels[nei[i]].belong) == 0) {
					std::vector<bool> bits = voxels[idx].hash.getBits();
					std::vector<bool> bits2 = voxels[idx].hash2.getBits();
					std::vector<bool> bits_ = voxels[nei[i]].hash.getBits();
					std::vector<bool> bits2_ = voxels[nei[i]].hash2.getBits();
					for (int b = 0; b < topo.knifes.size(); b++) {
						if (bits2[b] && bits2_[b] && bits[b] != bits_[b]) {
							Vector3 norm = topo.knifes[b].normal;
							bound.push_back(bits_[b] ? norm : -norm);
						}
					}
				}
				if (voxels[idx].itx > voxels[nei[i]].itx)bound.push_back(Vector3(-1, 0, 0));
				if (voxels[idx].itx < voxels[nei[i]].itx)bound.push_back(Vector3(1, 0, 0));
				if (voxels[idx].ity > voxels[nei[i]].ity)bound.push_back(Vector3(0, -1, 0));
				if (voxels[idx].ity < voxels[nei[i]].ity)bound.push_back(Vector3(0, 1, 0));
				if (voxels[idx].itz > voxels[nei[i]].itz)bound.push_back(Vector3(0, 0, -1));
				if (voxels[idx].itz < voxels[nei[i]].itz)bound.push_back(Vector3(0, 0, 1));
			}
		}
	}
	return bound;
}

std::vector<Vector3> Utility::genBound(std::vector<int> contact) {
	int lim = colthre;
	std::vector<Vector3> bound;
	if (contact[0] > lim)bound.push_back(Vector3(-1, 0, 0));
	if (contact[1] > lim)bound.push_back(Vector3(1, 0, 0));
	if (contact[2] > lim)bound.push_back(Vector3(0, -1, 0));
	if (contact[3] > lim)bound.push_back(Vector3(0, 1, 0));
	if (contact[4] > lim)bound.push_back(Vector3(0, 0, -1));
	if (contact[5] > lim)bound.push_back(Vector3(0, 0, 1));
	return bound;
}

Vector3 Utility::calGroupCenter(Group & group) {
	Vector3 cent(0,0,0);
	int num = 0;
	for (int i = 0; i < group.pieces.size(); i++) {
		for (int j = 0; j < pieces[group.pieces[i]].voxelsi.size(); j++) {
			int idx = pieces[group.pieces[i]].voxelsi[j];
			cent += voxels[idx].pos;
			num++;
		}
	}
	return cent / num;
}

float Utility::calGroupRate(Group & group) {
	Vector3 max(FLT_MIN, FLT_MIN, FLT_MIN);
	Vector3 min(FLT_MAX, FLT_MAX, FLT_MAX);
	for (int i = 0; i < group.pieces.size(); i++) {
		for (int j = 0; j < pieces[group.pieces[i]].voxelsi.size(); j++) {
			int idx = pieces[group.pieces[i]].voxelsi[j];
			max.x = std::max(max.x, voxels[idx].pos.x);
			max.y = std::max(max.y, voxels[idx].pos.y);
			max.z = std::max(max.z, voxels[idx].pos.z);
			min.x = std::min(min.x, voxels[idx].pos.x);
			min.y = std::min(min.y, voxels[idx].pos.y);
			min.z = std::min(min.z, voxels[idx].pos.z);
		}
	}
	Vector3 dis = (max - min);
	float mindis = std::min(dis.x, std::min(dis.y, dis.z));
	float maxdis = std::max(dis.x, std::max(dis.y, dis.z));
	return maxdis / mindis;
}

float Utility::max(std::vector<float> list) {
	float max = FLT_MIN;
	for (auto val : list) {
		max = std::max(max, val);
	}
	return max;
}

Worth Utility::calWorth(Group& group1, Group& group2) {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> rand(0, 1);
	//***************************************************
	/*
	Group group1_ = group1;
	Group group2_ = group2;
	topo.boundAva(group1_, genBound(calContact(group1_, group2)));
	topo.boundAva(group2_, genBound(calContact(group2_, group1)));
	int change = std::max((group1.avanum - group1_.avanum), (group2.avanum - group2_.avanum));
	*/
	//******************************************
	std::vector<int> contact = calContact(group1, group2);
	int contactface = 0;
	float contactarea = 0;
	float maxface = 0;
	for (auto c : contact) {
		maxface = std::max((float)c, maxface);
		contactarea += c;
		if (c > 0)contactface++;
	}
	float contactWarp = (contactarea / 6) / maxface;
	//*********************************************
	Group group = group1;
	/*
	group.borde = std::vector<int>(group1.bordenum + group2.bordenum);
	for (int i = 0; i < group1.bordenum; i++)group.borde[i] = group1.borde[i];
	for (int i = 0; i < group2.bordenum; i++)group.borde[i + group1.bordenum] = group2.borde[i];
	group.bordenum = group.borde.size();
	*/
	//*********************************************
	for (int i = 0; i<group2.pieces.size(); i++)appendPiece(group, pieces[group2.pieces[i]]);
	//checkCav(group, group2);
	int oriavanum = group.avanum;
	std::set<int> ids; ids.insert(group1.id); ids.insert(group2.id);
	//topo.boundAva(group, genBound(calContact(group, ids)));
	//******************************************
	appendGroup(group, group2);
	    
	int maxcoldis = INT_MIN;
	if (colmode && eachcolmode) {
		std::vector<int> coldis = calFarContact(group, ids);
		for (auto dis : coldis)maxcoldis = std::max(maxcoldis, dis);
	}
	//******************************************
	std::vector<std::set<int>> contactgroup;
	topo.boundAva(group, genBound(calBothContact(group1, group2, contactgroup)));
	checkCav(group);
	bool cav = false;
	cav = group.zcav;
	//******************************************
	float iscube = 0;
	int cubevolume = stx*sty*stz;
	if ((int)group1.volume%cubevolume == 0) iscube++;
	if ((int)group2.volume%cubevolume == 0) iscube++;
	float ctrl = 1;
	if (contact[4] > contlim || contact[5] > contlim) ctrl = 0;
	//******************************************
	/*
	std::vector<int> contact2 = calContact(group, ids);
	std::vector<int> contact1 = calBothContact(group1, group2);
	bool valid = false;
	for (int i = 0; i < 6; i++)if (contact2[i] != contact1[i])valid = true;
	if (valid) {
		for (int i = 0; i < 6; i++)std::cout << contact2[i] << " "; std::cout << std::endl;
		for (int i = 0; i < 6; i++)std::cout << contact1[i] << " "; std::cout << std::endl; std::cout << std::endl;
	}
	*/
	//****************************************************************
	//topo.boundFarAva(group, genBound(calFarContact(group, ids)));//has been moded
	//topo.boundAva(group, genBound(group, ids)); 
	//****************************************************************
	/*
	Hash hash; for (int i = 0; i < 6; i++)hash.addHash(true);
	for (int i = 0; i < group.pieces.size(); i++) {
		hash = hash & pieces[group.pieces[i]].hash;
	}
	std::vector<Vector3> bound;
	if (hash.getBit(0))bound.push_back(Vector3(-1, 0, 0));
	if (hash.getBit(1))bound.push_back(Vector3(1, 0, 0));
	if (hash.getBit(2))bound.push_back(Vector3(0, -1, 0));
	if (hash.getBit(3))bound.push_back(Vector3(0, 1, 0));
	if (hash.getBit(4))bound.push_back(Vector3(0, 0, -1));
	if (hash.getBit(5))bound.push_back(Vector3(0, 0, 1));
	topo.boundAva(group, bound);
	*/
	//******************************************
	//float addon = group1.trigger || group2.trigger ? 10000 : 0;
	//if (group.avanum == 0)return 0;
	//******************************************

	//******************************************
	int maxvolume = nx*ny*nz;
	int maxarea = 10000;
	//int volume = std::max(group1.volume, group2.volume);
	//int volume = std::min(group1.volume, group2.volume);
	int volume = group1.volume + group2.volume;
	std::vector<float> val;
	float lx = group.max.x - group.min.x;
	float ly = group.max.y - group.min.y;
	float lz = group.max.z - group.min.z;
	float rates[6] = { lx / ly, ly / lz, lz / lx , ly / lx, lz / ly,lx / lz };
	float maxrate = max(std::vector<float>(rates, rates + 6));
	//printf("%f\n", maxrate);
	float maxl = std::max(group.max.x - group.min.x, std::max(group.max.y - group.min.y, group.max.z - group.min.z));
	float minl = std::min(group.max.x - group.min.x, std::min(group.max.y - group.min.y, group.max.z - group.min.z));
	float rate = (float)group1.volume / group2.volume;
	if (rate > 1)rate = 1 / rate;
	//******************************************
	if (group.avanum == 0 || cav /* || contactarea < contlim|| maxrate > 5 */) { val.push_back(0); val.push_back(0); val.push_back(0); val.push_back(0); }
	else {
		val.push_back(group.avanum);
		//val.push_back(-iscube);
		val.push_back(contactface);
		//val.push_back(ctrl);
		val.push_back(maxarea - volume);
		//val.push_back(rate);
		//val.push_back(-lz);
		val.push_back(contactarea);
		//val.push_back(-maxrate);
		val.push_back(maxcoldis);//piority?????
	}
	//bound problem need to solve!!(完全平行於正交角時)
	return Worth(val);
	//******************************************
	//speed cav shape
	//speed ava

	//******************************************
	/*
	int incress = group.volume;
	incress -= group1.volume < group2.volume ? group1.volume : group2.volume;
	//******************************************
	int avadif = (group.avanum - group1.avanum) + (group.avanum - group2.avanum);
	//******************************************
	std::vector<int> contactfaces = calContact(group1, group2);
	int contactfacenum = 0;
	int contactnum = 0;
	for (int i = 0; i < 6; i++) {
		if (contactfaces[i] > 0)contactfacenum++;
		contactnum += contactfaces[i];
	}
	*/
	//******************************************
	/*
	float a1 = 0;
	float a2 = 0;
	float a3 = 0;
	float a4 = 1;
	float a5 = 0;
	return a1 * contactfacenum + a2 * contactnum +  a3 * incress + a4 * group.avanum + a5 * avadif;
	*/
}

Worth Utility::calWorth(Group& group1) {
	Group group = group1;
	std::set<int> ids; ids.insert(group1.id);
	//topo.boundAva(group, genBound(calContact(group, ids)));
	topo.boundAva(group, genBound(calContact(group)));
	
	//topo.boundFarAva(group, genBound(calFarContact(group, ids)));//has been moded
	//topo.boundAva(group, genBound(group, ids));
	std::vector<float> val;
	val.push_back(group.avanum);
	return Worth(val);
}

std::vector<int> Utility::calContact(Group & group1, Group & group2) {
	std::vector<int> contact = std::vector<int>(6, 0);
	if(group1.contactmap.count(group2.id)!=0)contact = group1.contactmap[group2.id];
	return contact;

	/*
	std::vector<int> contact2 = contact;
	contact = std::vector<int>(6, 0);
	for (auto p1 : group1.pieces) {
		for (auto p2 : group2.pieces) {
			for (int i = 0; i < 6; i++) {
				contact[i] += piececontact[p1][p2][i];
			}
		}
	}

	bool valid = false;
	for (int i = 0; i < 6; i++)if (contact2[i] != contact[i])valid = true;
	if (valid) {
		for (int i = 0; i < 6; i++)std::cout << contact2[i] << " "; std::cout << std::endl;
		for (int i = 0; i < 6; i++)std::cout << contact[i] << " "; std::cout << std::endl; std::cout << std::endl;
	}
	*/

	/*
	std::vector<int> belong(nx*ny*nz, 0);
	for (int i = 0; i < group1.bordenum; i++) {
		int idx = group1.borde[i];
		belong[idx] = 1;
	}
	for (int i = 0; i < group2.bordenum; i++) {
		int idx = group2.borde[i];
		belong[idx] = 2;
	}
	std::vector<int> contact= std::vector<int>(6, 0);
	for (int i = 0; i < group1.bordenum; i++) {
		int idx = group1.borde[i];
		int ix, iy, iz;
		voxels[idx].getXYZ(ix, iy, iz);
		int nei[6] = { -1, -1, -1, -1, -1, -1 };
		if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy)* nz + (iz));
		if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
		if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1) * nz + (iz));
		if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1) * nz + (iz));
		if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
		if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
		for (int i = 0; i < 6; i++) {
			if (nei[i] == -1)continue;
			if (belong[nei[i]] == 0)continue;
			if (belong[idx] != belong[nei[i]]) {
				contact[i]++;
			}
		}
	}
	return contact;
	*/

	/*
	for (int i = 0; i < group1.pieces.size();i ++) {
		for (int j = 0; j < pieces[group1.pieces[i]].voxelsi.size(); j++) {
			int idx = pieces[group1.pieces[i]].voxelsi[j];
			belong[idx] = 1;
		}
	}
	for (int i = 0; i < group2.pieces.size(); i++) {
		for (int j = 0; j < pieces[group2.pieces[i]].voxelsi.size(); j++) {
			int idx = pieces[group2.pieces[i]].voxelsi[j];
			belong[idx] = 2;
		}
	}
	for (int k = 0; k < group1.pieces.size(); k++) {
		for (int j = 0; j < pieces[group1.pieces[k]].voxelsi.size(); j++) {
			int idx = pieces[group1.pieces[k]].voxelsi[j];
			int ix, iy, iz;
			voxels[idx].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy)* nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1) * nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1) * nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
			for (int i = 0; i < 6; i++) {
				if (nei[i] == -1)continue;
				if (belong[nei[i]]==0)continue;
				if (belong[idx] != belong[nei[i]]) {
					contact[i]++;
				}
			}
		}
	}
	return contact;
	*/
}

std::vector<int> Utility::calContact(Group & group) {
	std::vector<int> contact(6, 0);
	for (auto iter : group.contactmap) {
		if (!idexist[iter.first])continue;
		if(groups[iter.first].removed)continue;
		if (iter.first == group.id)continue;
		for (int i = 0; i < 6; i++) {
			contact[i] += iter.second[i];
		}
	}
	return contact;
}

std::vector<int> Utility::calBothContact(Group & group1, Group & group2, std::vector<std::set<int>> & contactgroup) {
	std::vector<int> contact(6, 0);
	contactgroup = std::vector<std::set<int>>(6, std::set<int>());
	Group* grouplist[2] = { &group1, &group2 };
	for (int n = 0; n < 2; n++) {
		for (auto iter : grouplist[n]->contactmap) {
			if (!idexist[iter.first])continue;
			if(groups[iter.first].removed)continue;
			if (iter.first == group1.id || iter.first == group2.id)continue;
			for (int i = 0; i < 6; i++) {
				contact[i] += iter.second[i];
				if (iter.second[i] > 0)contactgroup[i].insert(iter.first);
			}
		}
	}
	return contact;
}

std::vector<int> Utility::calContact(Group & group, std::set<int> ids) {
	std::vector<int> contact(6, 0);
	for (int i = 0; i < group.pieces.size(); i++) {
		int cur = group.pieces[i];
		for (auto tar : pieces[cur].neighbor) {
			if (ids.count(pieces[tar].belong)==0 && !pieces[tar].removed) {
				for (int j = 0; j < 6; j++) {
					contact[j] += pieces[cur].contactmap[tar][j];
				}
			}
		}
	}
	//for (int i = 0; i < 6; i++)std::cout << contact[i] << " "; std::cout << std::endl;
	return contact;
	/*
	contact = std::vector<int>(6, 0);
	for (int i = 0; i < group.bordenum; i++) {
		int idx = group.borde[i];
		int ix, iy, iz;
		voxels[idx].getXYZ(ix, iy, iz);
		int nei[6] = { -1, -1, -1, -1, -1, -1 };
		if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy)* nz + (iz));
		if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
		if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1) * nz + (iz));
		if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1) * nz + (iz));
		if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
		if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
		for (int i = 0; i < 6; i++) {
			if (nei[i] == -1)continue;
			if (!voxels[nei[i]].exist)continue;
			if (voxels[nei[i]].removed)continue;
			if (ids.count(voxels[nei[i]].belong) == 0) {
				contact[i]++;
			}
		}
	}
	for (int i = 0; i < 6; i++)std::cout << contact[i] << " "; std::cout << std::endl; std::cout << std::endl;
	return contact;
	*/

	/*
	std::vector<int> contact(6, 0);
	for (int i = 0; i < 6; i++)contact[i] = 0;
	for (int k = 0; k < group.pieces.size(); k++) {
		for (int j = 0; j < pieces[group.pieces[k]].voxelsi.size(); j++) {
			int idx = pieces[group.pieces[k]].voxelsi[j];
			int ix, iy, iz;
			voxels[idx].getXYZ(ix, iy, iz);
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (ix - 1 >= 0)nei[0] = ((ix - 1) * ny * nz + (iy) * nz + (iz));
			if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy) * nz + (iz));
			if (iy - 1 >= 0)nei[2] = ((ix) * ny * nz + (iy - 1) * nz + (iz));
			if (iy + 1 < ny)nei[3] = ((ix) * ny * nz + (iy + 1) * nz + (iz));
			if (iz - 1 >= 0)nei[4] = ((ix) * ny * nz + (iy) * nz + (iz - 1));
			if (iz + 1 < nz)nei[5] = ((ix) * ny * nz + (iy) * nz + (iz + 1));
			for (int i = 0; i < 6; i++) {
				if (nei[i] == -1)continue;
				if (!voxels[nei[i]].exist)continue;
				if (voxels[nei[i]].removed)continue;
				if (ids.count(voxels[nei[i]].belong)==0) {
					contact[i]++;
				}
			}
		}
	}
	return contact;
	*/
}

std::vector<int> Utility::calFarContact(Group & group, std::set<int> ids) {
	int thre = 100;
	std::vector<int> contact(6, topo.maxcol);
	//for (int i = 0; i < 6; i++)if (group.avaset.count(i) == 0)contact[i] = 0;//only for simpmode
	int attemp = 2;
	if (attemp == 0) {
		for (int mode = 0; mode < 6; mode++) {
			if (group.avaset.count(mode) == 0) {contact[mode] = 0;continue;}
			int n, w1, w2;
			Vector3 dir;
			if (mode / 2 == 0) { n = nx; w1 = ny; w2 = nz; dir = Vector3(-1, 0, 0); }
			if (mode / 2 == 1) { n = ny; w1 = nx; w2 = nz; dir = Vector3(0, -1, 0); }
			if (mode / 2 == 2) { n = nz; w1 = nx; w2 = ny; dir = Vector3(0, 0, -1); }
			if (mode % 2 == 1) dir *= -1;
			std::vector<int> pixel(w1*w2); for (int j = 0; j < pixel.size(); j++) pixel[j] = -1;
			int mode_ = (mode / 2 * 2) + ((mode % 2 + 1) % 2);
			bool bk = false;
			std::vector<int> collist(n, 0);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < w1; j++) {
					for (int k = 0; k < w2; k++) {
						int idx = voxelId(i, j, k, mode_);
						bool same = ids.count(voxels[idx].belong) > 0;
						int dis = i - pixel[j * w2 + k];
						if (voxels[idx].exist && !voxels[idx].removed) {
							if (same)if(!voxels[idx].immo) pixel[j * w2 + k] = i;
							else if (pixel[j * w2 + k] >= 0) { collist[dis]++;}
						}
						if (pixel[j * w2 + k] >= 0 && !same && voxels[idx].immo) {
							collist[dis]+=10;
						}
					}
				}
			}
			for (int i = 0; i < n; i++) {
				if (collist[i] > thre) {
					contact[mode] = i;
					if (printdebug)printf("%d---%d ", i, collist[i]);
					break;
				}
			}
		}
		return contact;
	}
	else if(attemp == 1){
		std::set<int> idxts = group.idxts;
		//for (int i = 0; i < group.pieces.size(); i++) if (pieces[group.pieces[i]].volume > thre)idxts.insert(pieces[group.pieces[i]].it);
		for (auto idxt_cur : idxts) {
			int itx = idxt_cur / (nty * ntz);
			int ity = idxt_cur / ntz % nty;
			int itz = idxt_cur % ntz;
			for (int x = itx - 1; x >= 0; x--) {
				int dis = itx - x;
				if (dis >= contact[0])break;
				int idxt = x * nty * ntz + ity * ntz + itz;
				if (idxts.count(idxt) > 0)break;
				bool idxtbelong = false;
				for (auto p : supervoxels[idxt]) {
					if (pieces[p].volume <= thre)continue;
					int idx = pieces[p].voxelsi[0];
					if (voxels[idx].exist && !voxels[idx].removed && ids.count(voxels[idx].belong) == 0) {
						idxtbelong = true;
						contact[0] = dis;
						break;
					}
				}
				//if (!idxtbelong && supervoxeltouch[idxt])contact[0] = dis;
				if (!idxtbelong)contact[0] = dis;
				if (contact[0] > 0 && contact[0] < topo.maxcol)break;
			}
			for (int x = itx + 1; x < ntx; x++) {
				int dis = x - itx;
				if (dis >= contact[1])break;
				int idxt = x * nty * ntz + ity * ntz + itz;
				if (idxts.count(idxt) > 0)break;
				bool idxtbelong = false;
				for (auto p : supervoxels[idxt]) {
					if (pieces[p].volume <= thre)continue;
					int idx = pieces[p].voxelsi[0];
					if (voxels[idx].exist && !voxels[idx].removed && ids.count(voxels[idx].belong) == 0) {
						idxtbelong = true;
						contact[1] = dis;
						break;
					}
				}
				//if (!idxtbelong && supervoxeltouch[idxt])contact[1] = dis;
				if (!idxtbelong)contact[1] = dis;
				if (contact[1] > 0 && contact[1] < topo.maxcol)break;
			}
			for (int y = ity - 1; y >= 0; y--) {
				int dis = ity - y;
				if (dis >= contact[2])break;
				int idxt = itx * nty * ntz + y * ntz + itz;
				if (idxts.count(idxt) > 0)break;
				bool idxtbelong = false;
				for (auto p : supervoxels[idxt]) {
					if (pieces[p].volume <= thre)continue;
					int idx = pieces[p].voxelsi[0];
					if (voxels[idx].exist && !voxels[idx].removed && ids.count(voxels[idx].belong) == 0) {
						idxtbelong = true;
						contact[2] = dis;
						break;
					}
				}
				//if (!idxtbelong && supervoxeltouch[idxt])contact[2] = dis;
				if (!idxtbelong)contact[2] = dis;
				if (contact[2] > 0 && contact[2] < topo.maxcol)break;
			}
			for (int y = ity + 1; y < nty; y++) {
				int dis = y - ity;
				if (dis >= contact[3])break;
				int idxt = itx * nty * ntz + y * ntz + itz;
				if (idxts.count(idxt) > 0)break;
				bool idxtbelong = false;
				for (auto p : supervoxels[idxt]) {
					if (pieces[p].volume <= thre)continue;
					int idx = pieces[p].voxelsi[0];
					if (voxels[idx].exist && !voxels[idx].removed && ids.count(voxels[idx].belong) == 0) {
						idxtbelong = true;
						contact[3] = dis;
						break;
					}
				}
				//if (!idxtbelong && supervoxeltouch[idxt])contact[3] = dis;
				if (!idxtbelong)contact[3] = dis;
				if (contact[3] > 0 && contact[3] < topo.maxcol)break;
			}
			for (int z = itz - 1; z >= 0; z--) {
				int dis = itz - z;
				if (dis >= contact[4])break;
				int idxt = itx * nty * ntz + ity * ntz + z;
				if (idxts.count(idxt) > 0)break;
				bool idxtbelong = false;
				for (auto p : supervoxels[idxt]) {
					if (pieces[p].volume <= thre)continue;
					int idx = pieces[p].voxelsi[0];
					if (voxels[idx].exist && !voxels[idx].removed && ids.count(voxels[idx].belong) == 0) {
						idxtbelong = true;
						contact[4] = dis;
						break;
					}
				}
				//if (!idxtbelong && supervoxeltouch[idxt])contact[4] = dis;
				if (!idxtbelong)contact[4] = dis;
				if (contact[4] > 0 && contact[4] < topo.maxcol)break;
			}
			for (int z = itz + 1; z < ntz; z++) {
				int dis = z - itz;
				if (z - itz >= contact[5])break;
				int idxt = itx * nty * ntz + ity * ntz + z;
				bool idxtbelong = false;
				for (auto p : supervoxels[idxt]) {
					if (pieces[p].volume <= thre)continue;
					int idx = pieces[p].voxelsi[0];
					if (voxels[idx].exist && !voxels[idx].removed && ids.count(voxels[idx].belong) == 0) {
						idxtbelong = true;
						contact[5] = dis;
						break;
					}
				}
				//if (!idxtbelong && supervoxeltouch[idxt])contact[5] = dis;
				if (!idxtbelong)contact[5] = dis;
				if (contact[5] > 0 && contact[5] < topo.maxcol)break;
			}
		}
		return contact;
	}
	else {
		std::set<int> idxts;
		for (int i = 0; i < group.pieces.size(); i++) {
			if (pieces[group.pieces[i]].volume > thre)idxts.insert(pieces[group.pieces[i]].it);
		}
		for (int mode = 0; mode < 6; mode++) {
			if (group.avaset.count(mode) == 0) { contact[mode] = 0; continue; }
			int n, w1, w2, nt, wt1, wt2;
			Vector3 dir;
			if (mode / 2 == 0) { n = nx; w1 = ny; w2 = nz; nt = stx; wt1 = sty; wt2 = stz; dir = Vector3(-1, 0, 0); }
			if (mode / 2 == 1) { n = ny; w1 = nx; w2 = nz; nt = sty; wt1 = stx; wt2 = stz; dir = Vector3(0, -1, 0); }
			if (mode / 2 == 2) { n = nz; w1 = nx; w2 = ny; nt = stz; wt1 = stx; wt2 = sty; dir = Vector3(0, 0, -1); }
			if (mode % 2 == 1) dir *= -1;
			std::vector<int> pixel(w1*w2, -1);
			std::vector<bool> pixeldone(w1*w2, false);
			int mode_ = (mode / 2 * 2) + ((mode % 2 + 1) % 2);
			for (auto idxt_cur : idxts) {
				bool valid = true;
				//for (auto p : supervoxels[idxt_cur])if(ids.count(pieces[p].belong)==0)valid = false;
				if (valid && group.idxtcoldis[mode].count(idxt_cur) > 0) continue;
				std::vector<int> collist_(n, 0);
				int ti, tj, tk; voxelTId_(ti, tj, tk, mode_, idxt_cur);
				if (ti+1 < (n/nt) && idxts.count(voxelTId(ti+1, tj, tk, mode_)) > 0) continue;
				for (int i = ti * nt; i < n; i++) {
					int idxt = voxelTId(i / nt, tj, tk, mode_);
					
					bool jump = true;
					if (supervoxels[idxt].size() == 1) {
						if(idxt == idxt_cur){
							for (int j = tj * wt1; j < tj * wt1 + wt1; j++) for (int k = tk * wt2; k < tk * wt2 + wt2; k++) {
								pixel[j * w2 + k] = i + nt - 1;
							}
						}else{
							for (auto p : supervoxels[idxt]) {
								if (pieces[p].volume <= thre)continue;
								int idx = pieces[p].voxelsi[0];
								if (voxels[idx].exist && !voxels[idx].removed && ids.count(voxels[idx].belong) == 0) {
									jump = false;
								}
							}
						}
						
					}else jump = false;
					if (idxt != idxt_cur && jump) { i += nt - 1; continue; }
					
					int pixeldonecnt = 0;
					for (int j = tj * wt1; j < tj * wt1 + wt1; j++) {
						for (int k = tk * wt2; k < tk * wt2 + wt2; k++) {
							if (pixeldone[j * w2 + k])continue;
							int idx = voxelId(i, j, k, mode_);
							bool same = ids.count(voxels[idx].belong) > 0;
							int dis = i - pixel[j * w2 + k];
							if (voxels[idx].exist && !voxels[idx].removed) {
								if (same)if (!voxels[idx].immo) pixel[j * w2 + k] = i;
								else if (pixel[j * w2 + k] >= 0) { collist_[dis]++; pixeldone[j * w2 + k] = true; }
							}
							if (pixel[j * w2 + k] >= 0 && !same && voxels[idx].immo) { collist_[dis] += 10; pixeldone[j * w2 + k] = true;}
							if (pixeldone[j * w2 + k])pixeldonecnt++;
						}
						if (pixeldonecnt == wt1*wt2)break;
					}
					if (pixeldonecnt == wt1*wt2)break;
				}
				group.idxtcoldis[mode][idxt_cur] = collist_;
			}
			std::vector<int> collist(n, 0);
			for (auto pair : group.idxtcoldis[mode]) {
				for (int i = 0; i < n; i++) {
					collist[i] += pair.second[i];
				}
			}
			for (int i = 0; i < n; i++) {
				if (collist[i] > thre) {
					contact[mode] = i;
					break;
				}
			}
		}
		return contact;
	}
}

//********************************************************

void Utility::outputPiece() {
	for (int i = 0; i < pieces.size(); i++) {
		Mesh mesh = pieces[i].mesh;
		char str1[20];
		char str2[20];
		sprintf(str1, "piece_%d", i);
		iglMachine.put(str1, mesh.vertices, mesh.indices);
		sprintf(str2, "pool\\piece_%d.obj", i);
		iglMachine.writeFile(str1, str2);
		std::cout << "Write " << str1 << std::endl;
	}
}

void Utility::outputKnife() {
	FILE* fp = fopen("knifeinfo.txt", "w");
	for (int i = 0; i < topo.knifeIdx.size(); i++) {
		fprintf(fp, "%d\n", topo.knifeIdx[i]);
	}
	fclose(fp);
}

void Utility::outputZip() {
	//float fix = 1.0f;//????????????????????????????????????
	float fix = outputfix;
	std::vector<std::vector<int>> pieceingrid(ntx*nty*ntz);
	for (int i = 0; i < pieces.size(); i++) {
		if (pieces[i].volume == stx*sty*stz)pieces[i].iscube = true;
		pieceingrid[voxels[pieces[i].voxelsi[0]].it].push_back(i);
	}
	//*************************************
	FILE* fp = fopen("pieceinfo.txt", "w");
	fprintf(fp, "%d\n", groups_final.size());
	for (int i = 0; i < groups_final.size(); i++) {
		Group group = groups[groups_final[i]];
		//*************************************
		int minus = 0;
		std::vector<bool> dispear(group.pieces.size(), true);
		if (clearmode) {
			for (int j = 0; j < group.pieces.size(); j++) {
				Piece piece = pieces[group.pieces[j]];
				if (!piece.iscube)dispear[j] = false;
				int it = voxels[piece.voxelsi[0]].it;
				int itx, ity, itz; voxels[piece.voxelsi[0]].getTXYZ(itx, ity, itz);
				int nei[6] = { -1, -1, -1, -1, -1, -1 };
				if (itx - 1 >= 0)nei[0] = ((itx - 1) * nty * ntz + (ity)* ntz + (itz));
				if (itx + 1 < ntx)nei[1] = ((itx + 1) * nty * ntz + (ity)* ntz + (itz));
				if (ity - 1 >= 0)nei[2] = ((itx)* nty * ntz + (ity - 1) * ntz + (itz));
				if (ity + 1 < nty)nei[3] = ((itx)* nty * ntz + (ity + 1) * ntz + (itz));
				if (itz - 1 >= 0)nei[4] = ((itx)* nty * ntz + (ity)* ntz + (itz - 1));
				if (itz + 1 < ntz)nei[5] = ((itx)* nty * ntz + (ity)* ntz + (itz + 1));
				for (int k = 0; k < 6; k++) {
					if (nei[k] == -1)continue;
					if (pieceingrid[it].size()>1)continue;
					if (!pieces[pieceingrid[nei[k]][0]].iscube) {
						for (auto p : pieceingrid[nei[k]]) {
							if (pieces[p].belong == pieces[pieceingrid[it][0]].belong) {
								dispear[j] = false;
							}
						}
					}
				}
				if (itz == 0 || itz == ntz - 1)dispear[j] = false;
				if (dispear[j])minus++;
			}
		}
		else {
			dispear = std::vector<bool>(group.pieces.size(), false);
		}
		//*************************************
		fprintf(fp, "%d\n", group.pieces.size()- minus);
		for (int j = 0; j < group.pieces.size(); j++) {
			if (dispear[j])continue;
			Piece piece = pieces[group.pieces[j]];
			int itx, ity, itz; voxels[piece.voxelsi[0]].getTXYZ(itx, ity, itz);
			Vector3 pos = ld + Vector3(itx * stx * l, ity * sty * l, itz * stz * l) + Vector3(stx * l / 2, sty * l / 2, stz * l / 2);
			Vector3 size_ = Vector3(stx * l / 2, sty * l / 2, stz * l / 2);
			Vector3 ld_ = pos - size_, ru_ = pos + size_;
			Vector3 ld__ = ld_ + Vector3(fix, fix, fix), ru__ = ru_ - Vector3(fix, fix, fix);
			ld_ = ld__; ru_ = ru__;
			int nei[6] = { -1, -1, -1, -1, -1, -1 };
			if (itx - 1 >= 0)nei[0] =  ((itx - 1) * nty * ntz + (ity) * ntz + (itz));
			if (itx + 1 < ntx)nei[1] = ((itx + 1) * nty * ntz + (ity) * ntz + (itz));
			if (ity - 1 >= 0)nei[2] =  ((itx) * nty * ntz + (ity - 1) * ntz + (itz));
			if (ity + 1 < nty)nei[3] = ((itx) * nty * ntz + (ity + 1) * ntz + (itz));
			if (itz - 1 >= 0)nei[4] =  ((itx) * nty * ntz + (ity) * ntz + (itz - 1));
			if (itz + 1 < ntz)nei[5] = ((itx) * nty * ntz + (ity) * ntz + (itz + 1));
			for (int k = 0; k < topo.edgenum; k++) {
				for (auto p : group.pieces) {
					if (voxels[pieces[p].voxelsi[0]].it == nei[k]) {
						if (k == 0)ld_.x = ld__.x - fix;
						if (k == 1)ru_.x = ru__.x + fix;
						if (k == 2)ld_.y = ld__.y - fix;
						if (k == 3)ru_.y = ru__.y + fix;
						if (k == 4)ld_.z = ld__.z - fix;
						if (k == 5)ru_.z = ru__.z + fix;
					}
				}
			}
			pos = (ld_ + ru_) / 2;
			Vector3 size = (ru_ - ld_);
			//*******************************************
			if (!manhmode) {
				fprintf(fp, "%f %f %f %f %f %f\n", pos.x, pos.y, pos.z, size.x, size.y, size.z);
			}
			else {
				Cube manhcube = manhCubeMap[voxels[piece.voxelsi[0]].hash4];
				fprintf(fp, "%f %f %f %f %f %f\n", manhcube.cent.x, manhcube.cent.y, manhcube.cent.z, manhcube.size.x, manhcube.size.y, manhcube.size.z);
			}
			std::vector<bool> bits = voxels[piece.voxelsi[0]].hash.getBits();
			std::vector<bool> bits2 = voxels[piece.voxelsi[0]].hash2.getBits();
			for (int k = 0; k < topo.edgenum; k++) {
				if (bits2[k]) {
				//if (bits2[k] && topo.knifeExist[k]) {
					bool valid = true;
					for (int m = 0; m < group.pieces.size(); m++) {
						if (j == m)continue;
						Voxel voxel = voxels[pieces[group.pieces[m]].voxelsi[0]];
						if (voxel.hash2.getBit(k) && bits[k] != voxel.hash.getBit(k)) {
							valid = false;
							break;
						}
					}
					Vector3 cent = topo.knifes[k].center;
					Vector3 dir = topo.knifes[k].normal;
					if (!bits[k])dir *= -1;
					if (valid) {
						cent += dir.normalize() * fix;
					}
					fprintf(fp, "%f %f %f ", cent.x, cent.y, cent.z);
					fprintf(fp, "%f %f %f ", dir.x, dir.y, dir.z);
				}
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
}

void Utility::outputPiece_voxel() {
	int cnt = 0;
	for (int i = 0; i < pieces.size(); i++) {
		char filename[20]; sprintf(filename, "pool\\shape_%d.voxel", cnt);
		std::ofstream ofs = std::ofstream(filename, std::ios::binary);
		std::vector<Vector3> poss;
		for (int j = 0; j < pieces[i].voxels.size(); j++) {
			poss.push_back(pieces[i].voxels[j]);
		}
		int ps = poss.size();
		ofs.write((char*)&ps, sizeof(int));
		ofs.write((char*)&l, sizeof(float));
		for (auto pos : poss) {
			ofs.write((char*)&pos.x, sizeof(float));
			ofs.write((char*)&pos.y, sizeof(float));
			ofs.write((char*)&pos.z, sizeof(float));
		}
		ofs.close();
		if (pieces[i].mesh.indices.size() > 0) {
			char objname[20]; sprintf(objname, "shape_%d", cnt);
			outputMesh(objname, pieces[i].mesh);
		}
		std::cout << "Write " << "shape_" << cnt << ".voxel" << std::endl;
		cnt++;
	}
	FILE* fp = fopen("shapeinfo.txt", "w");
	fprintf(fp, "%d\n", groups.size());
	for (int i = 0; i < groups.size(); i++) {
		Vector3 dir = Vector3(0, 0, 0);
		if (groups[i].avanum > 0) {
			for (int j = 0; j < groups[i].avanum; j++) dir += groups[i].ava[j];
			dir /= groups[i].avanum;
			printf("%f %f %f\n", dir.x, dir.y, dir.z);
			fprintf(fp, "%f %f %f\n", dir.x, dir.y, dir.z);
		}
		else fprintf(fp, "0 0 0\n");
	}
	fclose(fp);
}

void Utility::outputGroup_voxel() {
	int cnt = 0;
	std::map<int, int> cntmap;
	for (int i = 0; i < groups_final.size(); i++) {
		cntmap[groups_final[i]] = cnt;
		Group group = groups[groups_final[i]];
		char filename[20]; sprintf(filename, "pool\\shape_%d.voxel", cnt);
		std::vector<Vector3> poss;
		for (int j = 0; j < group.pieces.size(); j++) {
			for (int k = 0; k < pieces[group.pieces[j]].voxels.size(); k++) {
				poss.push_back(pieces[group.pieces[j]].voxels[k]);
			}
		}
		int ps = poss.size();
		bool on = true;
		if (on) {
			std::ofstream ofs = std::ofstream(filename, std::ios::binary);
			ofs.write((char*)&ps, sizeof(int));
			ofs.write((char*)&l, sizeof(float));
			for (auto pos : poss) {
				ofs.write((char*)&pos.x, sizeof(float));
				ofs.write((char*)&pos.y, sizeof(float));
				ofs.write((char*)&pos.z, sizeof(float));
			}
			ofs.close();
		}
		std::cout << "Write " << "shape_" << cnt << ".voxel" << std::endl;
		cnt++;
	}
	FILE* fp = fopen("shapeinfo.txt", "w");
	fprintf(fp, "%d\n", cnt);
	for (int i = 0; i < groups_final.size(); i++) {
		Group group = groups[groups_final[i]];
		Vector3 dir = Vector3(0, 0, 0);
		if (group.avanum > 0) {
			
			//for (int j = 0; j < group.avanum; j++)dir += group.ava[j];
			//dir /= group.avanum;
			//fprintf(fp, "%f %f %f\n", dir.x, dir.y, dir.z);
			//int a = group.avanum / 2;
			//fprintf(fp, "%f %f %f\n", group.ava[a].x, group.ava[a].y, group.ava[a].z);
			
			topo.renewAva(group);
			if (colmode) {
				fprintf(fp, "%f %f %f %d\n", group.maxDir.x, group.maxDir.y, group.maxDir.z, group.maxDirDis);
				printf("%f %f %f : %d\n", group.maxDir.x, group.maxDir.y, group.maxDir.z, group.maxDirDis);
			}
			else {
				int n = group.ava.size() / 2;
				fprintf(fp, "%f %f %f 0\n", group.ava[n].x, group.ava[n].y, group.ava[n].z);
				printf("%f %f %f : 0\n", group.ava[n].x, group.ava[n].y, group.ava[n].z);
			}
		}
		else {
			fprintf(fp, "0 0 0 0\n");
			printf("0 0 0 0\n");
		}
	}
	for (int i = 0; i < groups_final_assem.size(); i++) {
		fprintf(fp, "%d\n", cntmap[groups_final_assem[i]]);
	}
	fclose(fp);
}

float Utility::outputEnergy() {
	int num = groups_final.size();
	//*****************************
	float maxlen = topo.maxcol;
	float coldis = 0;
	for (int i = 0; i < groups_final.size(); i++) {
		Group group = groups[groups_final[i]];
		topo.renewAva(group);
		if (group.maxDirDis >= maxlen || group.maxDirDis < 0)group.maxDirDis = maxlen;
		coldis += 1-((float)group.maxDirDis / maxlen);
	}
	coldis /= num;
	//*****************************
	std::vector<float> volumes;
	float sum = 0;
	float max = 0;
	for (int i = 0; i < groups_final.size(); i++) {
		Group group = groups[groups_final[i]];
		volumes.push_back(group.volume);
		max = std::max(max, group.volume);
		sum += group.volume;
	}
	float stddif = 0;
	float mean = (sum / max) / num;
	for (auto vol : volumes) {
		float rate = vol / max;
		stddif += (rate - mean)*(rate - mean);
	}
	stddif = std::sqrtf(stddif / num);
	//*****************************
	energy[0] = (float)num;
	energy[1] = coldis;
	energy[2] = stddif;
	//std::cout << num << " " << coldis << " " << stddif << std::endl;
	std::cout << energy[0] << " " << energy[1] << " " << energy[2] << std::endl;
	std::cout << groups[groups_final[groups_final.size() - 1]].volume / (nx*ny*nx) << std::endl;
	float E = energy[0] + 10 * energy[1] + energy[2];
	std::cout << E << std::endl;
	//17.1042
	//19.9273
	return E;
}

void Utility::outputGroup() {
	int cnt = 0;
	for (int i = 0; i < groups.size(); i++) {
		if (!idexist[groups[i].id])continue;
		Mesh mesh = groups[i].getMesh(pieces);//need ref
		//mesh.merge();
		char str1[20];
		char str2[20];
		sprintf(str1, "shape_%d", cnt);
		iglMachine.put(str1, mesh.vertices, mesh.indices);
		sprintf(str2, "pool\\shape_%d.obj", cnt);
		iglMachine.writeFile(str1, str2);
		std::cout << "Write " << str1 << std::endl;
		cnt++;
	}
	FILE* fp = fopen("shapeinfo.txt","w");
	fprintf(fp, "%d\n", cnt);
	for (int i = 0; i < groups.size();i++) {
		if (!idexist[groups[i].id])continue;
		Vector3 dir = Vector3(0, 0, 0);
		if (groups[i].avanum > 0) {
			for (int j = 0; j < groups[i].avanum; j++)dir += groups[i].ava[j];
			dir /= groups[i].avanum;
			fprintf(fp, "%f %f %f\n", dir.x, dir.y, dir.z);
		}
		else fprintf(fp, "0 0 0\n");
	}
	fclose(fp);
}

void Utility::outputMesh(char * str, Mesh mesh) {
	char path[50]; sprintf(path, "pool\\%s.obj", str);
	FILE* fp = fopen(path, "w");
	for (int i = 0; i < mesh.vertices.size(); i+=3) {
		fprintf(fp, "v %f %f %f\n", mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]);
	}
	for (int i = 0; i < mesh.indices.size(); i += 3) {
		fprintf(fp, "f %d %d %d\n", mesh.indices[i] + 1, mesh.indices[i + 1] + 1, mesh.indices[i + 2] + 1);
	}
	fclose(fp);
}

//parallel problem?