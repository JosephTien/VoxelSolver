#ifndef MESH_H
#define MESH_H
#include <pch.h>
#include <geometry.h>
//****************************
/*hide upon lite
#include "cgaltool.h"
*/
//****************************
//replace upon lite
class CGALTool {
public:
	CGALTool() {}
};
//****************************
class Mesh
{
public:

	Mesh() {}
    Mesh(std::vector<float> vertices, std::vector<unsigned int> indices){
        this->vertices = vertices;
        this->indices = indices;
        this->normals = vertices;//quick assign
    }
    std::vector<float> vertices;
    std::vector<unsigned int> indices;
	std::vector<float> normals;
	std::set<int> pick;
	std::vector<std::vector<unsigned int>> neighbor;
	CGALTool cgaltool;
    //******************************************************************
	float dotProduct(Vector3 a, Vector3 b);
    Vector3 getVertice(int idx);
    Vector3 getNormal(int idx);
    unsigned int pushVertice(Vector3 v);
    void putVertice(int idx, Vector3 v);
	void pushIndice(unsigned int a, unsigned int b, unsigned int c);
	void pushIndice(unsigned int a, unsigned int b, unsigned int c, unsigned int d);
	void purne();
	void rotateTo(Vector3 vec);
    Mesh slice(Vector3 c,Vector3 n);
    void regenByPlateIntersec(Vector3 c, Vector3 n);
    Mesh genCube(Vector3 ld, Vector3 ru);
    Mesh genTri(Vector3 c, float r);
	void filte(Vector3 c, Vector3 n);
	void addMesh(Mesh mesh);
	void merge();
	void removeFace(Vector3 n);
	static float area(Vector3 a, Vector3 b, Vector3 c);

	//******************************************************************
	//replace upon lite
	void convexHull() {};
	//******************************************************************
	/*hide upon lite
	void load(const char* filename);
	void load_rc(QFile * qfile);
	void loadObj(const char* filename, std::vector<float> &verts, std::vector<float> &norms, std::vector<unsigned int> &facs);
	void loadObj_rc(QFile * qfile, std::vector<float> &verts, std::vector<float> &norms, std::vector<unsigned int> &facs);
	void fillHole();
	void convexHull();
	void simplify();
	*/
	
    //******************************************************************
    //******************************************************************
	/*Bkup
	 *
    Plane cuttingPlane;
    std::vector<std::vector<unsigned int>> neighbor;
    std::set<Edge> edges;
    std::vector<int> detourIdxs;
    std::vector<int> detourIdxs_all;
    int detourSPIdx;
    //******************************************************************
	void fillHole(Vector3 n);
    int calDetourByPlane(Vector3 c,Vector3 n);
    void cutByDetour(int state);
    bool Mesh::checkDetour();
    void calNeighbor();
    void fix();
    void regenNormals();
    Vector3 detourCenter();
    void sortDetour();

	void calNeighbor();
	*/

};

#endif // MESH_H
