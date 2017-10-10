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
	std::vector<Piece> pieces;
	std::vector<Group> groups;
	std::vector<int> groups_final;
	std::priority_queue<GroupLink> groupLink;
	std::vector<int> groupIdxMap;//std::map<int, int> groupIdxMap;
	std::vector<bool> idexist;
	int nx, ny, nz;
	int groupidcnt = 0;
	Vector3 ld, ru;
	const float l = 2.0f;
	float far = 100.0f;// 10.0f;
	clock_t t1, t2;
	//*************************
	void tic() {
		t1 = clock();
	}
	void toc() {
		t2 = clock(); printf("%lf sec\n", (t2 - t1) / (double)(CLOCKS_PER_SEC)); t1 = clock();
	}
	//*************************
	void genRandomTest(int k);
	void genPiece(std::string filename);
	void genPiece(std::string filename, bool output);
	void genPieceGroupMesh(std::string filename);
	//*************************
	void calBound();
	void genPiece();
	void genPiece_voxel();
	void appendPiece(Group& group, Piece& piece);
	void MergeGroup(Group& group1, Group& group2);
	void initGroup();
	void calTouch(int p);
	void initLink();
	void initLink_voxel();
	float calWorth(Group& group1, Group& group2);
	std::vector<int> calContact(Group & group, std::set<int> ids);
	int calContact(Group & group1, Group & group2);
	void optimize();
	void iterate();
	//*************************
	int voxelId(int i, int j, int k, int mode);
	void genPiece_voxel_bfs();
	void genVoxelByKnife();
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
	void genVoxelOutput();
	void collectLast();
	//*************************
	void outputPiece();
	void outputPiece_voxel();
	void outputGroup();
	void outputGroup_voxel();
	void outputMesh(char * str, Mesh mesh);
	void outputKnife();
	
};

#endif // UTILITY_H
