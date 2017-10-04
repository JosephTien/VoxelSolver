#include "utility.h"

float gr() {
	return (((float)rand() / RAND_MAX) - 0.5f) * 2;
}

void Utility::genPieceGroupMesh(std::string filename){
	bool realMerge = false;
    std::ifstream ifs;
    ifs.open(filename);
    std::string line;
    Vector3 ld, ru; float x,y,z; int n,l,m;
    char str[20];
    std::stringstream ss;
    std::getline( ifs, line);ss = std::stringstream(line);
    ss >> x >> y >> z; ld = Vector3(x,y,z);
    ss >> x >> y >> z; ru = Vector3(x,y,z);
    std::cout << "left low : " << ld << std::endl;
    std::cout << "right up : " << ru << std::endl;
    //********************************************************
    std::vector<Vector3> cs,ns;
    std::getline( ifs, line);ss = std::stringstream(line);
    ss >> n >> l;
    std::cout << "knife number : " << n << std::endl;
    std::cout << "shape number : " << l << std::endl;
    for(int j=0;j<n;j++){
        std::getline( ifs, line);ss = std::stringstream(line);
        ss >> x >> y >> z;cs.push_back(Vector3(x,y,z));
        ss >> x >> y >> z;ns.push_back(Vector3(x,y,z));
        std::cout << "knife_" << j << " : " << cs[j] << " " << ns[j] <<std::endl;
    }
    //********************************************************
    std::cout << "________________________________________" << std::endl;
    iglMachine.reset();
    iglMachine.command("NEW null");
    for(int j=0;j<n;j++){
        //Mesh mesh = Mesh().genCube(cs[j]-Vector3(100000,100000,0), cs[j]+Vector3(100000,100000,100000));
        Mesh mesh = Mesh().genTri(cs[j], 100000);
        mesh.rotateTo(ns[j]); //notice the reverse case
        sprintf(str,"knife_%d", j);
        iglMachine.put(str, mesh.vertices, mesh.indices);
    }
    //********************************************************
    int cnt=0;
    for(int j=0;j<l;j++){
        std::cout << "________________________________________" << std::endl;
        std::getline( ifs, line);ss = std::stringstream(line);
        ss >> m;
        int tar = cnt++;
        std::cout << "Generate shape_" << tar << " with " << m << " pieces\n>" << std::endl;
        sprintf(str,"NEW shape_%d", tar);iglMachine.command(str);
        for(int k=0;k<m;k++){
            std::getline( ifs, line);
            Mesh mesh = Mesh().genCube(ld, ru);
            iglMachine.put("temp", mesh.vertices, mesh.indices);
            const char * chs = line.c_str();
            int len = strlen(chs);
            for(int i=0;i<len;i++){
                sprintf(str,"%c temp temp knife_%d", chs[i], i);iglMachine.command(str);
            }
            sprintf(str,"ADD shape_%d temp", tar);iglMachine.command(str);
        }
		if (realMerge) { sprintf(str, "+ shape_%d shape_%d null", tar, tar); iglMachine.command(str); }
        sprintf(str,"WRITE shape_%d", tar);iglMachine.command(str);
    }
    std::cout << "________________________________________\nDone!" << std::endl;
	//********************************************************
}

void Utility::genRandomTest(int k) {
	FILE * fp = fopen("randomTest.txt", "w");
	fprintf(fp, "-1 -1 -1 1 1 1\n");
	fprintf(fp, "%d 0\n", k);
	for (int i = 0; i < k; i++) {
		fprintf(fp, "%f %f %f %f %f %f\n", gr(), gr(), gr(), gr(), gr(), gr());
	}
	fclose(fp);
}

int Utility::voxelId(int i, int j, int k,int mode) {
	int idx;
	if (mode == 0)idx = i * ny * nz + j * nz + k;
	if (mode == 1)idx = (nx-i-1) * ny * nz + j * nz + k;
	if (mode == 2)idx = j * ny * nz + i * nz + k;
	if (mode == 3)idx = j * ny * nz + (ny-i-1) * nz + k;
	if (mode == 4)idx = j * ny * nz + k * nz + i;
	if (mode == 5)idx = j * ny * nz + k * nz + (nz-i-1);
	return idx;
}

void Utility::voxelDirSearch(int mode) {
	Piece piece;
	int p = pieces.size();
	int n, w1, w2;
	//*******************************************
	if (mode/2==0) { n = nx; w1 = ny; w2 = nz; }
	if (mode/2==1) { n = ny; w1 = nx; w2 = nz; }
	if (mode/2==2) { n = nz; w1 = nx; w2 = ny; };

	//*******************************************
	std::vector<int> pixel(w1*w2);for (int j = 0; j < w1; j++) for (int k = 0; k < w2; k++) pixel[j*w2 + k] = 0;
	bool nfirst = false;
	for (int i = 0; i < n; i++) {
		nfirst = true;
		int anything = false;
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
		if (anything==0 && !nfirst)break;
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

//*******************************************

void Utility::voxelBfs() { //將所有touchedge先放入，並一組一組分別長
	std::set<int> vs;
	for (int i = 0; i < voxels.size(); i++) vs.insert(i);
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
void Utility::genVoxelSeen() {//將moveable direction塞入hash
	clock_t t1, t2; t1 = clock();
	//*****************************
	std::vector<Capsule> caps = std::vector<Capsule>();
	for (int i = 0; i < topo.edgenum; i++) {
		topo.knifes.push_back(Plane(topo.splitNorm[i], topo.getEdgeCent(i)));
		topo.knifeIdx.push_back(i);
		caps.push_back(Capsule(topo.vertices[topo.edges[i].ia], topo.vertices[topo.edges[i].ib], topo.radii));
	}
	std::vector<bool> bits = std::vector<bool>(topo.knifes.size() * 3);
	nx = (ru.x - ld.x) / l + 1;
	ny = (ru.y - ld.y) / l + 1;
	nz = (ru.z - ld.z) / l + 1;
	voxels = std::vector<Voxel>(nx*ny*nz);
	//*******************************************
	#pragma omp parallel for
	for (int i = 0; i < voxels.size(); i++) {
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
			if (caps[j].collide(pos, l / 2, between))immo = true;
			if (caps[j].collide(pos, 1.5 *l, between))voxels[idx].touchid.insert(j);
		}
		voxels[idx].setXYZ(ix, iy, iz);
		voxels[idx].exist = false;
		voxels[idx].immo = immo;
	}
	//*******************************************
	#pragma omp parallel for
	for (int i = 0; i < voxels.size(); i++) {
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
	for (int i = 0; i < voxels.size(); i++)hs.insert(voxels[i].hash);
	std::cout << "Label number : " << hs.size() << std::endl;
	//*******************************************
	t2 = clock(); printf("%lf sec\n", (t2 - t1) / (double)(CLOCKS_PER_SEC)); t1 = clock();
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

void Utility::voxelDirSeen() {
	voxelDirSeen(4);
	voxelDirSeen(5);
	voxelDirSeen(0);
	voxelDirSeen(1);
	voxelDirSeen(2);
	voxelDirSeen(3);
}

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

void Utility::voxelBfsSeen() {
	std::set<int> last;
	for (int i = 0; i < voxels.size(); i++)last.insert(i);

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

void Utility::genVoxelOutput() {
	Piece piece;
	for (int i = 0; i < voxels.size(); i++) {
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
	//********************************************
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
				if (voxels[nei[j]].belong != n)continue;
				samecnt++;
			}
			if (samecnt < 6) {
				pieces[n].voxels.push_back(voxels[i].pos);
			}
		}
	}
}

void Utility::voxelCollectSeen(int mode) {//group只用來記錄方向
	Piece piece;
	for (int i = 0; i < voxels.size(); i++) {
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
	clock_t t1, t2; t1 = clock();

	std::vector<Capsule> caps = std::vector<Capsule>();
	for (int i = 0; i < topo.edgenum; i++) {
		topo.knifes.push_back(Plane(topo.splitNorm[i], topo.getEdgeCent(i)));
		topo.knifeIdx.push_back(i);
		caps.push_back(Capsule(topo.vertices[topo.edges[i].ia], topo.vertices[topo.edges[i].ib], topo.radii));
	}
	std::vector<bool> bits = std::vector<bool>(topo.knifes.size() * 3);
	nx = (ru.x - ld.x) / l + 1;
	ny = (ru.y - ld.y) / l + 1;
	nz = (ru.z - ld.z) / l + 1;
	voxels = std::vector<Voxel>(nx*ny*nz);
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
					bits[i * 3] = caps[i].collide(pos, l, between);//it is L
					bits[i * 3 + 1] = between;
					bits[i * 3 + 2] = (pos - topo.knifes[i].center).dot(topo.knifes[i].normal) >= 0;
					if (bits[i * 3]) voxels[idx].touchid.insert(i);
					if (caps[i].collide(pos, l/2, between))immo = true;
				}
				voxels[idx].hash.assign(bits);
				voxels[idx].setXYZ(ix, iy, iz);
				voxels[idx].exist = false;
				voxels[idx].immo = immo;
			}
		}
	}
	t2 = clock(); printf("%lf sec\n", (t2 - t1) / (double)(CLOCKS_PER_SEC)); t1 = clock();
	//**************************************
	//test
	/* 
	pieces = std::vector<Piece>(topo.edgenum);
	for (int i = 0; i < voxels.size(); i++) {
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
		std::cout << pieces[i].voxels.size() << std::endl;
	}
	*/
}

void Utility::genVoxelByKnife() {
	nx = (ru.x - ld.x) / l + 1;
	ny = (ru.y - ld.y) / l + 1;
	nz = (ru.z - ld.z) / l + 1;
	voxels = std::vector<Voxel>(nx*ny*nz);
	#pragma omp parallel for
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
				int idx = ix * ny * nz + iy * nz + iz;
				std::vector<bool> bits = std::vector<bool>(topo.knifes.size());
				for (int i = 0; i < topo.knifes.size(); i++) {
					bits[i] = topo.knifes[i].distanceToPoint(pos) >= 0;
				}
				voxels[idx] = Voxel(Hash().assign(bits), pos, idx);
				voxels[idx].setXYZ(ix, iy, iz);
			}
		}
	}
	//************************************************************************************
	#pragma omp parallel for
	for (int i = 0; i < voxels.size(); i++) {
		int ix = i / (ny * nz);
		int iy = i / nz % ny;
		int iz = i % nz;
		int nei[6] = { -1, -1, -1, -1, -1, -1 };
		if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
		if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
		if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
		if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
		if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
		if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
		int samecnt = 0;
		for (int j = 0; j < 6; j++) {
			if (nei[j] != -1 && voxels[nei[j]].hash == voxels[i].hash) {
				samecnt++;
			}
		}
		if (samecnt == 6)voxels[i].exist = false;
	}
}

void Utility::genPiece(std::string filename) {
	genPiece(filename, true);
}

void Utility::genPiece(std::string filename, bool output) {
	std::ifstream ifs;
	ifs.open(filename);
	std::string line;
	Vector3 ld, ru; float x, y, z; int n, l, m;
	char str[20];
	std::stringstream ss;
	std::getline(ifs, line); ss = std::stringstream(line);
	ss >> x >> y >> z; ld = Vector3(x, y, z);
	ss >> x >> y >> z; ru = Vector3(x, y, z);
	std::cout << "left low : " << ld << std::endl;
	std::cout << "right up : " << ru << std::endl;
	//********************************************************
	std::vector<Vector3> cs, ns;
	std::getline(ifs, line); ss = std::stringstream(line);
	ss >> n >> l;
	std::cout << "knife number : " << n << std::endl;
	std::cout << "shape number : " << l << std::endl;
	for (int j = 0; j<n; j++) {
		std::getline(ifs, line); ss = std::stringstream(line);
		ss >> x >> y >> z; cs.push_back(Vector3(x, y, z)); 
		ss >> x >> y >> z; ns.push_back(Vector3(x, y, z));
		topo.knifes.push_back(Plane(ns[ns.size() - 1], cs[cs.size() - 1]));
		std::cout << "knife_" << j << " : " << cs[j] << " " << ns[j] << std::endl;

	}
	//********************************************************
	std::cout << "                    " << std::endl;
	pieces.push_back(Piece(Mesh().genCube(ld, ru)));
	for (int j = 0; j< n; j++) {
		std::cout << "                    " << std::endl;
		std::cout << "Apply Knife " << j+1 << "/" << n << std::endl;
		std::vector<Piece> pieces_next = std::vector<Piece>();
		#pragma omp parallel for ordered
		for (int i = 0; i < pieces.size(); i++) {
			Piece piece = pieces[i];
			Piece piece_ = Piece(piece.mesh.slice(cs[j], ns[j]), piece.hash);
			#pragma omp ordered
			if(true){
				if (piece.mesh.vertices.size() > 0) {
					piece.hash.addHash(true);
					pieces_next.push_back(piece);
				}
				if (piece_.mesh.vertices.size() > 0) {
					piece_.hash.addHash(false);
					pieces_next.push_back(piece_);
				}
			}
		}
		pieces.swap(pieces_next);
		std::cout << "Pieces Num : " << pieces.size() << std::endl;
		std::cout << "____________________" << std::endl;
	}
	//*****************************
	/* merge test
	while (pieces.size() > 1) {
		std::vector<Piece> pieces_next;
		for (int i = 0; i < pieces.size(); i+=2) {
			if (i + 1 >= pieces.size())pieces_next.push_back(pieces[i]);
			else {
				pieces[i].mesh.addMesh(pieces[i + 1].mesh);
				pieces_next.push_back(pieces[i]);
			}
		}
		pieces.swap(pieces_next);
		std::cout << "Pieces Num : " << pieces.size() << std::endl;
	}
	*/
	//*****************************
	if (!output)return;
	for (int i = 0; i < pieces.size(); i++) {
		Mesh mesh = pieces[i].mesh;
		char str1[20];
		char str2[20];
		sprintf(str1, "shape_%d", i);
		iglMachine.put(str1, mesh.vertices, mesh.indices);
		sprintf(str2, "pool\\shape_%d.obj", i);
		iglMachine.writeFile(str1, str2);
		std::cout << "Write " << str1 << std::endl;
	}
}

void Utility::genPiece_voxel() {
	clock_t t1, t2; t1 = clock();
	//*****************************
	genVoxelByKnife();
	t2 = clock(); printf("%lf sec\n", (t2 - t1) / (double)(CLOCKS_PER_SEC)); t1 = clock();
	//************************************************************************************
	std::map<Hash, Piece> hashMap;
	for (int i = 0; i < voxels.size(); i++) {
		int ix = i / (ny * nz);
		int iy = i / nz % ny;
		int iz = i % nz;
		Vector3 pos = ld + Vector3(l*ix, l*iy, l*iz);
		if (hashMap.count(voxels[i].hash) == 0) {
			hashMap[voxels[i].hash] = Piece(Mesh(), voxels[i].hash);//insert
		}
		if (voxels[i].exist)hashMap[voxels[i].hash].voxels.push_back(pos);
		hashMap[voxels[i].hash].voxelsi.push_back(i);
		hashMap[voxels[i].hash].volume++;
	}
	std::map<Hash, Piece>::iterator it;
	for (it = hashMap.begin(); it != hashMap.end(); it++) {
		pieces.push_back(it->second);
	}
	t2 = clock(); printf("%lf sec\n", (t2 - t1) / (double)(CLOCKS_PER_SEC)); t1 = clock();
	//************************************************************************************
	#pragma omp parallel for
	for (int i = 0; i < pieces.size(); i++) {
		//voxel mode gen .voxel but no mesh
		//for (int j = 0; j < pieces[i].voxels.size(); j++) {
		//	Vector3 pos = pieces[i].voxels[j];
		//	pieces[i].mesh.addMesh(Mesh().genCube(pos - Vector3(l / 2, l / 2, l / 2), pos + Vector3(l / 2, l / 2, l / 2)));//it is L
		//}
		for (int j = 0; j < pieces[i].voxelsi.size(); j++) {
			voxels[pieces[i].voxelsi[j]].belong = i;
		}
		pieces[i].id = i;
	}
	t2 = clock(); printf("%lf sec\n", (t2 - t1) / (double)(CLOCKS_PER_SEC)); t1 = clock();
	//************************************************************************************
	//printf("l = %f, size = %d x %d x %d, cost %lf sec\n", l, nx, ny, nz, (t2 - t1) / (double)(CLOCKS_PER_SEC));
	//************************************************************************************
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

void Utility::initGroup() {
	groups = std::vector<Group>(pieces.size());
	#pragma omp parallel for
	for (int i = 0; i < pieces.size(); i++) {
		Group group;
		group.initAva();
		topo.calTouch(pieces[i]);
		appendPiece(group, pieces[i]);
		group.id = i;
		groups[i] = group;
	}
	groupidcnt = pieces.size();
	groupIdxMap = std::vector<int>(groups.size());
	for (int i = 0; i < groups.size(); i++) {
		groupIdxMap[groups[i].id] = i;//groupIdxMap.insert(std::pair<int, int>(groups[i].id, i));
		idexist.push_back(true);
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
				groupLink.push(GroupLink(groups[i].id, groups[j].id, calWorth(groups[i], groups[j]), groups[i].volume + groups[j].volume));
			}
		}
	}
}

void Utility::initLink_voxel() {
	#pragma omp parallel for
	for (int i = 0; i < voxels.size(); i++) {
		int ix = i / (ny * nz);
		int iy = i / nz % ny;
		int iz = i % nz;
		int nei[6] = { -1, -1, -1, -1, -1, -1 };
		if (ix - 1 >= 0)nei[0] = ((ix - 1)* ny * nz + (iy)* nz + (iz));
		if (ix + 1 < nx)nei[1] = ((ix + 1) * ny * nz + (iy)* nz + (iz));
		if (iy - 1 >= 0)nei[2] = ((ix)* ny * nz + (iy - 1)* nz + (iz));
		if (iy + 1 < ny)nei[3] = ((ix)* ny * nz + (iy + 1)* nz + (iz));
		if (iz - 1 >= 0)nei[4] = ((ix)* ny * nz + (iy)* nz + (iz - 1));
		if (iz + 1 < nz)nei[5] = ((ix)* ny * nz + (iy)* nz + (iz + 1));
		int samecnt = 0;
		for (int j = 0; j < 6; j++) {//assume piece id is equal to group id
			if (nei[j] != -1) {
				int a = voxels[i].belong;
				int b = voxels[nei[j]].belong;
				groups[a].neighbor.insert(groups[b].id);
				groups[b].neighbor.insert(groups[a].id);
			}
		}
	}
	for (int i = 0; i < groups.size();i++) {
		for (auto jj : groups[i].neighbor) {
			int j = groupIdxMap[jj];
			if (i >= j)continue;
			groupLink.push(GroupLink(groups[i].id, groups[j].id, calWorth(groups[i], groups[j]), groups[i].volume + groups[j].volume));
		}
	}
}

void Utility::MergeGroup(Group& group1, Group& group2) {
	Group group = group1;
	for (int i = 0; i < group2.pieces.size(); i++) {
		appendPiece(group, pieces[group2.pieces[i]]);
	}
	group.id = groupidcnt++;
	idexist.push_back(true);
	idexist[group1.id] = idexist[group2.id] = false;
	group.neighbor.swap(std::set<int>());
	for (auto id : group1.neighbor)if (id != group1.id && id != group2.id && idexist[id]) { group.neighbor.insert(id);}
	for (auto id : group2.neighbor)if (id != group1.id && id != group2.id && idexist[id]) { group.neighbor.insert(id);}
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
	//*******************************
	
	groups.push_back(group);
	groupIdxMap.push_back((int)(groups.size() - 1));//groupIdxMap.insert(std::pair<int, int>(group.id, (int)(groups.size() - 1)));
	for (auto id : group.neighbor) {
		groups[groupIdxMap[id]].neighbor.insert(group.id);
		float worth = calWorth(group, groups[groupIdxMap[id]]);
		float volume = (group.volume + groups[groupIdxMap[id]].volume);
		groupLink.push(GroupLink(group.id, id, worth, volume));
	}
	
}

void Utility::appendPiece(Group& group, Piece& piece){
	group.pieces.push_back(piece.id);
	group.volume += piece.volume;
	topo.purneAva(group, piece);
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
	ld -= Vector3(radii * 2, radii * 2, radii * 2);
	ru += Vector3(radii * 2, radii * 2, radii * 2);
	std::cout << ld << std::endl;
	std::cout << ru << std::endl;
}

void Utility::optimize() {
	int cnt = groups.size();
	while (true) {
		GroupLink gl = groupLink.top();
		groupLink.pop();
		//std::cout << "Top link energy = " << gl.worth << std::endl;
		if (!idexist[gl.ida]|| !idexist[gl.idb])continue; //if (groupIdxMap.count(gl.ida) == 0 || groupIdxMap.count(gl.idb) == 0)continue;
		if (gl.worth == 0)break;
		MergeGroup(groups[groupIdxMap[gl.ida]], groups[groupIdxMap[gl.idb]]);
		cnt--;
		std::cout << "Group number left : " << cnt << std::endl;
		if (cnt <= 1)break;
	}
}

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
	for (int i = 0; i < groups.size(); i++) {
		if (!idexist[groups[i].id])continue;
		char filename[20]; sprintf(filename, "pool\\shape_%d.voxel", cnt);

		std::ofstream ofs = std::ofstream(filename, std::ios::binary);
		std::vector<Vector3> poss;
		for (int j = 0; j < groups[i].pieces.size(); j++) {
			for (int k = 0; k < pieces[groups[i].pieces[j]].voxels.size(); k++) {
				poss.push_back(pieces[groups[i].pieces[j]].voxels[k]);
			}
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
		std::cout << "Write " << "shape_" << cnt << ".voxel" << std::endl;
		cnt++;
	}
	FILE* fp = fopen("shapeinfo.txt", "w");
	fprintf(fp, "%d\n", cnt);
	for (int i = 0; i < groups.size(); i++) {
		if (!idexist[groups[i].id])continue;
		Vector3 dir = Vector3(0, 0, 0);
		if (groups[i].avanum > 0) {
			for (int j = 0; j < groups[i].avanum; j++)dir += groups[i].ava[j];
			dir /= groups[i].avanum;
			printf("%f %f %f\n", dir.x, dir.y, dir.z);
		}
		else fprintf(fp, "0 0 0\n");
	}
	fclose(fp);
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

float Utility::calWorth(Group& group1, Group& group2) {
	//return 1;
	Group group = group1;
	for (int i = 0; i < group2.pieces.size(); i++) {
		topo.purneAva(group, pieces[group2.pieces[i]]);
	}
	return group.avanum;
}
//********************************************************

//旋轉至最少knife
//旋轉改最大值

/* BKUP : csg method
void Utility::genPiece(std::string filename) {
std::ifstream ifs;
ifs.open(filename);
std::string line;
Vector3 ld, ru; float x, y, z; int n, l, m;
char str[20];
std::stringstream ss;
std::getline(ifs, line); ss = std::stringstream(line);
ss >> x >> y >> z; ld = Vector3(x, y, z);
ss >> x >> y >> z; ru = Vector3(x, y, z);
std::cout << "left low : " << ld << std::endl;
std::cout << "right up : " << ru << std::endl;
//********************************************************
std::vector<Vector3> cs, ns;
std::getline(ifs, line); ss = std::stringstream(line);
ss >> n >> l;
std::cout << "knife number : " << n << std::endl;
std::cout << "shape number : " << l << std::endl;
for (int j = 0; j<n; j++) {
std::getline(ifs, line); ss = std::stringstream(line);
ss >> x >> y >> z; cs.push_back(Vector3(x, y, z));
ss >> x >> y >> z; ns.push_back(Vector3(x, y, z));
std::cout << "knife_" << j << " : " << cs[j] << " " << ns[j] << std::endl;
}
//********************************************************
std::cout << "________________________________________" << std::endl;
cgalnef.genCube(ld, ru);
std::cout << "Cube generated..." << std::endl;
cgalnef.genPiece(cs, ns);
std::cout << "Pieces generated..." << std::endl;
for (int i = 0; i < cgalnef.Ns.size(); i++) {
Mesh mesh;
cgaltool.readFromOFFStream(mesh.vertices, mesh.indices, cgalnef.getOffStr(cgalnef.Ns[i]));
char str1[20], str2[20];
sprintf(str1, "shape_%d", i);
iglMachine.put(str1, mesh.vertices, mesh.indices);
sprintf(str2, "pool\\shape_%d.obj", i);
iglMachine.writeFile(str1, str2);
std::cout << "Pieces " << i << " Writed!" << std::endl;
}
}
*/

/* BKUP : my slicing step
 *
if(chs[i]=='*')mesh.slice(cs[i], ns[i]);
if(chs[i]=='-')mesh.slice(cs[i], ns[i]*-1);
fillHole(mesh);
*/

/* BKUP : dividing slicing
 *
std::vector<Mesh> meshes;
std::vector<Mesh> meshes_next;
meshes.push_back(Mesh().genCube(ld, ru));
for(int j=0;j<n;j++){
ss = std::stringstream(line);
ss >> x >> y >> z;
Vector3 c(x,y,z);
ss >> x >> y >> z;
Vector3 n(x,y,z);
for(int i=0;i<meshes.size();i++){
Mesh mesh = meshes[i];
Mesh mesh_ = mesh.slice(c, n);
fillHole(mesh);fillHole(mesh_);
if(mesh.indices.size()>0)meshes_next.push_back(mesh);
if(mesh_.indices.size()>0)meshes_next.push_back(mesh_);
}
meshes.swap(meshes_next);
meshes_next.swap(std::vector<Mesh>());
}
int cnt=0;
for(int i=0;i<meshes.size();i++){
Mesh mesh = meshes[i];
char str1[20];
char str2[20];
sprintf(str1, "shape_%d", i);
iglMachine.put(str1, mesh.vertices, mesh.indices);
sprintf(str2, "pool\\shape_%d.obj", i);
iglMachine.writeFile(str1,str2);
}
*/
