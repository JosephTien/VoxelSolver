#include "datastructure.h"

void Group::initAva() {
	int num = 0;
	int n = 360 / (int)plus;
	ava = std::vector<Vector3>(n*n*n);
	avalist = std::vector<Vector3>(n*n*n);
	int cnt = 0;
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
	avanum = ava.size();
}

void Group::initAvaLev() {
	int s = 360 / plus;
	avalevel = std::vector<int>(s*s*s, 0);
}

Mesh Group::getMesh(std::vector<Piece> & piecesref) {
	Mesh mesh;
	for (int i = 0; i < pieces.size(); i++) {
		mesh.addMesh(piecesref[pieces[i]].mesh);
	}
	return mesh;
}