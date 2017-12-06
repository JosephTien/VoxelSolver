#ifndef UTILITY_H
#define UTILITY_H

#include "pch.h"
#include "topo.h"
#include "time.h"
//****************************************************************
/*hide upon lite
#include "iglmachine.h"
#include "cgalnef.h"
*/
//****************************************************************
//replace upon lite
class CGALNef {
public:
	CGALNef() {}
};
class IglMachine {
public:
	IglMachine() {}
	void put(std::string, std::vector<float>, std::vector<unsigned int>) {}
	void writeFile(std::string, std::string) {}
	void command(std::string) {}
	void reset() {}
};
//****************************************************************

class Utility
{
public:
    Utility(){}
	CGALNef cgalnef;
	IglMachine iglMachine;
	Topo topo;
	std::vector<Voxel> voxels;
	std::vector<std::vector<int>> tnear;
	std::vector<std::vector<bool>> tbits2;
	std::vector<std::vector<int>> supervoxels;//
	std::vector<bool> supervoxeltouch;//
	std::vector<std::vector<int>> supervoxelblock;
	std::map<int, std::map<int, std::vector<int>>> piececontact;
	
	std::vector<Piece> pieces;
	std::vector<Group> groups;
	std::set<int> finalgroupset;
	std::vector<int> groups_final;
	std::vector<int> groups_final_assem;
	std::priority_queue<GroupLink> groupLink;
	std::vector<int> groupIdxMap;//std::map<int, int> groupIdxMap;
	std::vector<bool> idexist;
	int nx, ny, nz;
	int groupidcnt = 0;
	Vector3 ld, ru;
	float energy[3];

	const float l = 2.0f;
	float far = 10.0f;// 10.0f;
	int levlim = 3;
	int colthre = 10;
	int contlim = 10;
	bool tagmode = false;
	bool nearmode = false;
	bool manhmode = false;
	bool clearmode = false;
	bool checkcavmode = true;//似乎有一些問題，顯示voxel看看?
	bool quickcheckcav = true;
	bool colmode = true;
	bool eachcolmode = false;
	bool printdebug = true;
	bool printdebug_l2 = false;
	
	std::map<Hash, Cube> manhCubeMap;
	float manh = topo.radii + topo.radii;
	Vector3 move;
	int stx = 12, sty = 12, stz = 12;
	int ntx, nty, ntz;
	clock_t t1, t2;
	//*************************
	void tic() {
		t1 = clock();
	}
	void toc() {
		t2 = clock(); printf("%lf sec\n", (t2 - t1) / (double)(CLOCKS_PER_SEC));
	}
	//*************************
	void initVar() {
		voxels = std::vector<Voxel>();
		supervoxels = std::vector<std::vector<int>>();
		supervoxeltouch = std::vector<bool>();
		supervoxelblock = std::vector<std::vector<int>>();
		piececontact = std::map<int, std::map<int, std::vector<int>>>();
		pieces = std::vector<Piece>();
		groups = std::vector<Group>();
		finalgroupset = std::set<int>();
		groups_final = std::vector<int>();
		groups_final_assem = std::vector<int>();
		groupLink = std::priority_queue<GroupLink>();
		groupIdxMap = std::vector<int>();
		idexist = std::vector<bool>();
		groupidcnt = 0;
	}
	void genRandomTest(int k);
	void genPiece(std::string filename);
	void genPiece(std::string filename, bool output);
	void genPieceGroupMesh(std::string filename);
	//*************************
	void calBound();
	void genPiece();
	void genPiece_voxel();
	void appendPiece(Group& group, Piece& piece);
	void appendGroup(Group& group, Group& group2);
	void collectCav(Group& group);
	void checkCav(Group& group);
	void checkCav(Group& group, Piece& piece);
	void checkCav(Group& group, Group& group2);
	void MergeGroup(Group& group1, Group& group2);
	void initGroup();
	void initGroupBorde(Group & group);
	void purneGroupBorde(Group & group);
	void calTouch(int p);
	std::vector<Vector3> genBound(std::vector<int> contact);
	void initLink();
	void initLink_voxel();
	Worth calWorth(Group& group1, Group& group2);
	Worth calWorth(Group& group1);
	std::vector<Vector3> genBound(Group & group, std::set<int> ids);
	Vector3 calGroupCenter(Group & group);
	float calGroupRate(Group & group);
	std::vector<int> calContact(Group & group, std::set<int> ids);
	std::vector<int> calFarContact(Group & group, std::set<int> ids);
	std::vector<int> calContact(Group & group1, Group & group2);
	std::vector<int> calContact(Group & group);
	std::vector<int> calBothContact(Group & group1, Group & group2, std::vector<std::set<int>> & contactgroup);
	void optimize();
	void removeGroup(int tar);
	void iterate();
	void recalAssem();
	void noopt();
	void purneAvaByNei(int i);
	//*************************
	int voxelId(int i, int j, int k, int mode);
	int voxelTId(int i, int j, int k, int mode);
	void voxelTId_(int & i, int & j, int & k, int mode, int idxt);
	void genSuperVoxelBlock();
	void genVoxelByKnife_autotune();
	float genVoxelByKnife();
	void preview();
	void previewVoxelByKnife(std::vector<Plane> & knifes);
	void previewPiece_voxel(std::vector<Plane> & knifes);
	void genPiece_voxel_bfs();
	void genVoxel();
	void genVoxelSeen();
	void calVoxelSeen(int mode);
	void voxelDirSeen();
	void voxelDirSeen(int mode);
	void voxelDirSearch(int mode);
	bool voxelBfs(int p, std::set<int>& vs);
	void voxelBfs(int p);
	void voxelBfs();
	void voxelBfsOnce(int p, std::vector<int>& queue, int& flag, std::set<int>& touchedge, int mode);
	void voxelBfsSeen();
	bool voxelBfsSeenOnce(std::set<int> & last, int cmpbit);
	//*************************
	void voxelCollectSeen(int mode);
	void genVoxelMesh(int p, int mode);
	void purneVoxel(Group & group);
	void genVoxelOutput();
	void collectLast();
	//*************************
	void outputPiece();
	void outputZip();
	void outputPiece_voxel();
	void outputGroup();
	void outputGroup_voxel();
	float outputEnergy();
	void outputMesh(char * str, Mesh mesh);
	void outputKnife();
	//*************************
	float max(std::vector<float> list);
};

#endif // UTILITY_H
