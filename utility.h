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
	std::priority_queue<GroupLink> groupLink;
	std::vector<int> groupIdxMap;//std::map<int, int> groupIdxMap;
	std::vector<bool> idexist;
	int nx, ny, nz;
	int groupidcnt = 0;
	Vector3 ld, ru;
	const float l = 2.0f;
	void genRandomTest(int k);
    void genPieceGroupMesh(std::string filename);
	void calBound();
	void genVoxel();
	void genVoxelSeen();
	void calVoxelSeen(int state);
	void voxelDirSeen();
	void voxelDirSeen(int state);
	void voxelBfsSeen();
	bool voxelBfsSeenOnce(std::set<int> & last, int cmpbit);
	void voxelCollectSeen(int state);
	void genVoxelMesh(int p, int state);
	void genVoxelByKnife();
	void genPiece_voxel();
	void genPiece();
	void genPiece(std::string filename);
	void genPiece(std::string filename, bool output);
	void appendPiece(Group& group, Piece& piece);
	void MergeGroup(Group& group1, Group& group2);
	void initGroup();
	void initLink();
	void initLink_voxel();
	void optimize();
	int voxelId(int i, int j, int k, int state);
	void voxelDirSearch(int state);
	void voxelBfsOnce(int p, std::vector<int>& queue, int& flag, std::set<int>& touchedge, int state);
	bool voxelBfs(int p, std::set<int>& vs);
	void voxelBfs(int p);
	void voxelBfs();
	void outputMesh(char str[50], Mesh mesh);
	void outputPiece();
	void genVoxelOutput();
	void outputPiece_voxel();
	void outputGroup();
	void outputGroup_voxel();
	void outputKnife();
	float calWorth(Group& group1, Group& group2);
};

#endif // UTILITY_H
