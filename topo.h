#ifndef TOPO_H
#define TOPO_H
#include "pch.h"
#include "geometry.h"
#include "datastructure.h"

class TopoEdge {
public:
	TopoEdge(int a, int b) {
		ia = a;
		ib = b;
	}
	int ia;
	int ib;
	float fixa;
	float fixb;
};

class AngArg{
public:
	AngArg() {
		val = 0;
	}
	AngArg(std::vector<float> angles, float val) {
		this->angles.swap(angles);
		this->val = val;
	}
	std::vector<float> angles;
	float val;
	int val2;
	bool operator < (const AngArg& tar) const
	{
		return val < tar.val;
	}
	bool operator == (const AngArg& tar) const
	{
		return val > tar.val;
	}
	bool operator > (const AngArg& tar) const
	{
		return val > tar.val;
	}
	std::vector<float> difto(std::vector<float> tar) {
		for (int i = 0; i < tar.size(); i++) {
			tar[i] = tar[i] - angles[i];
		}
		return tar;
	}
	void adddif(std::vector<float> tar, float mul) {
		for (int i = 0; i < tar.size(); i++) {
			angles[i] = angles[i] + tar[i] * mul;
			if (angles[i] > 179)angles[i] = 179;
			if (angles[i] < 0)angles[i] = 0;
		}
	}
	void adddif(std::vector<float> tar) {
		adddif(tar, 1);
	}
};

class Topo
{
public:
	Topo() {
	
	}
	std::vector<Vector3> vertices;
	std::vector<TopoEdge> edges;
	std::vector<Vector3> splitNorm;
	std::vector<Vector3> splitNorm_ori;
	std::vector<Plane> knifes;
	std::vector<int> knifeIdx;//用edge idx 代替表達output的knife
	std::vector<std::set<int>> verticeedge;
	std::vector<std::set<int>> verticevertice;
	std::vector<bool> isimpo;
	std::vector<float> angles;
	std::vector<Capsule> caps;
	std::vector<Vector3> avalist;
	int optMode = 0;
	const float radii = 5.0f;
	int verticenum;
	int edgenum;
	int lesstype = 0;
	void genKnife();
	void genAllKnife();
	void genCapsule();
	std::vector<Plane> genKnife(std::vector<Vector3> splitNorm, std::vector<int>& knifeIdx);
	void read();
	float angleFix(float angle);
	void fixAngleDistance();
	Vector3 getEdgeCent(int idx) { return (vertices[edges[idx].ia] + vertices[edges[idx].ib]) / 2; }
	Vector3 getEdgeVec(int idx) { return (vertices[edges[idx].ib] - vertices[edges[idx].ia]).normalize(); }
	void prepareData();
	float calAngleArgVal(std::vector<float> angles);
	void calTouch(Piece & piece);
	void calTouchBound(Piece & piece);
	void purneAva(Group &group, Piece &piece);
	void boundAva(Group &group, std::vector<Vector3>);
	void boundAva_angles(Group &group, std::vector<Vector3>);
	void boundFarAva(Group &group, std::vector<Vector3>);
	std::vector<float> genRandomAngles();
	std::vector<float> genRandomAngles(float from, float to, float lev);
	void outputRotateArg();
	Vector3 genRandomOrgdir();
	float Topo::calPerAngle(Vector3 norm, Vector3 vec, Vector3 dir);
	void beeOpt();
	void geneOpt();
	void geneOpt(float from, float to, float lev);
	void geneOpt(std::vector<bool> pos, float from, float to, float lev);
	void geneOpt_org();
	void geneLess();
	void geneOrg();
	void atomOpt();
	Vector3 calNorm(int i, float angle);
	float calCrossVal(int i, Vector3 norm, float disthre);
	std::vector<Vector3> calTouchBound(std::set<TouchInfo> &tis);
};

#endif // TOPO_H