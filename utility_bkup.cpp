#include "utility.h"

float gr() {
	return (((float)rand() / RAND_MAX) - 0.5f) * 2;
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

void Utility::genPieceGroupMesh(std::string filename) {
	bool realMerge = false;
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
	iglMachine.reset();
	iglMachine.command("NEW null");
	for (int j = 0; j<n; j++) {
		//Mesh mesh = Mesh().genCube(cs[j]-Vector3(100000,100000,0), cs[j]+Vector3(100000,100000,100000));
		Mesh mesh = Mesh().genTri(cs[j], 100000);
		mesh.rotateTo(ns[j]); //notice the reverse case
		sprintf(str, "knife_%d", j);
		iglMachine.put(str, mesh.vertices, mesh.indices);
	}
	//********************************************************
	int cnt = 0;
	for (int j = 0; j<l; j++) {
		std::cout << "________________________________________" << std::endl;
		std::getline(ifs, line); ss = std::stringstream(line);
		ss >> m;
		int tar = cnt++;
		std::cout << "Generate shape_" << tar << " with " << m << " pieces\n>" << std::endl;
		sprintf(str, "NEW shape_%d", tar); iglMachine.command(str);
		for (int k = 0; k<m; k++) {
			std::getline(ifs, line);
			Mesh mesh = Mesh().genCube(ld, ru);
			iglMachine.put("temp", mesh.vertices, mesh.indices);
			const char * chs = line.c_str();
			int len = strlen(chs);
			for (int i = 0; i<len; i++) {
				sprintf(str, "%c temp temp knife_%d", chs[i], i); iglMachine.command(str);
			}
			sprintf(str, "ADD shape_%d temp", tar); iglMachine.command(str);
		}
		if (realMerge) { sprintf(str, "+ shape_%d shape_%d null", tar, tar); iglMachine.command(str); }
		sprintf(str, "WRITE shape_%d", tar); iglMachine.command(str);
	}
	std::cout << "________________________________________\nDone!" << std::endl;
	//********************************************************
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
		std::cout << "Apply Knife " << j + 1 << "/" << n << std::endl;
		std::vector<Piece> pieces_next = std::vector<Piece>();
#pragma omp parallel for ordered
		for (int i = 0; i < pieces.size(); i++) {
			Piece piece = pieces[i];
			Piece piece_ = Piece(piece.mesh.slice(cs[j], ns[j]), piece.hash);
#pragma omp ordered
			if (true) {
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
