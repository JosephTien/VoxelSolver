#include "datastructure.h"

void Group::initAva() {
	for (float x = 0; x < 360; x += plus) {
		for (float y = 0; y < 360; y += plus) {
			for (float z = 0; z < 360; z += plus) {
				Matrix4 mat;
				mat.rotateX(x); mat.rotateY(y); mat.rotateZ(z);
				ava.push_back(mat*Vector3(1, 0, 0));
			}
		}
	}
	avanum = ava.size();
}

Mesh Group::getMesh(std::vector<Piece> & piecesref) {
	Mesh mesh;
	for (int i = 0; i < pieces.size(); i++) {
		mesh.addMesh(piecesref[pieces[i]].mesh);
	}
	return mesh;
}