#include "datastructure.h"

void Group::initAva() {
	int num = 0;
	if (!simpmode) {
		int n = 360 / (int)plus;
		avalist = std::vector<Vector3>(n*n*n);
		ava = std::vector<Vector3>(n*n*n);
		for (float x = 0.1; x < 360; x += plus) {
			for (float y = 0.1; y < 360; y += plus) {
				for (float z = 0.1; z < 360; z += plus) {
					Matrix4 mat;
					mat.rotateX(x); mat.rotateY(y); mat.rotateZ(z);
					ava[num] = mat*Vector3(1, 0, 0);
					avalist[num] = ava[num];
					avaset.insert(num);
					num++;
				}
			}
		}
	}
	else {
		avalist = std::vector<Vector3>(6);
		ava = std::vector<Vector3>(6);
		Vector3 samp[6] = { Vector3(-1,0,0) ,Vector3(1,0,0), Vector3(0,-1,0), Vector3(0,1,0), Vector3(0,0,-1),Vector3(0,0,1) };
		for (int i = 0; i < 6; i++) {
			//samp[i] -= Vector3(0.001f, 0.001f, 0.001f);
			ava[num] = samp[i];
			avalist[num] = ava[num];
			avaset.insert(num);
			num++;
		}
	}
	avanum = ava.size();
}

void Group::initAvaLev(int maxcol) {
	if (!simpmode) {
		int s = 360 / plus;
		avalevel = std::vector<int>(s*s*s, maxcol);
	}
	else {
		avalevel = std::vector<int>(6, maxcol);
	}
	coldis = std::vector<int>(6, maxcol);
}

Mesh Group::getMesh(std::vector<Piece> & piecesref) {
	Mesh mesh;
	for (int i = 0; i < pieces.size(); i++) {
		mesh.addMesh(piecesref[pieces[i]].mesh);
	}
	return mesh;
}