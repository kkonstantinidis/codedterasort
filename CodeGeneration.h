#ifndef _CMR_CODEGENERATION
#define _CMR_CODEGENERATION

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <unordered_map>
#include <assert.h>

using namespace std;

typedef vector< int > Block;
typedef vector< Block > ParallelClass;
typedef set< int > InputSet;
typedef int BlockId;
typedef vector< BlockId > ParallelClassIds;
typedef set< int > NodeSet;
typedef unsigned int SGId;
typedef pair< int, int > Vpair;//< destId, inputId >
typedef vector< Vpair > VpairList;

class CodeGeneration {
private:
	int N;
	int K;
	int k;
	int q;
	map< int, InputSet > M;//ranks of user will increase from block to block and then from parallel class to parallel class continuously
	vector< ParallelClassIds > resolutionIds;
	vector< NodeSet > shufflingGroups;
	map< NodeSet, SGId > SGIdMap;
	map< NodeSet, int > commonSubfile;
	map< int, vector< NodeSet > > NodeSGMap;
	map< int, NodeSet > FileNodes;
	
public:
	CodeGeneration(int _N, int _K, int _k, int _q);
	~CodeGeneration(){}
	InputSet& getM(int nodeId){return M[nodeId];}
	map< NodeSet, SGId >& getSGIdMap(){return SGIdMap;}
	NodeSet& getNodeSetFromFileID(unsigned long fid){return FileNodes[fid];}
	vector< NodeSet >& getNodeSGContain(int nid){return NodeSGMap[nid];}
	SGId getSGId(NodeSet ns);
	int getCommonSubfileId(NodeSet ns){return commonSubfile[ns];}
	static void printNodeSet(NodeSet ns);

private:
	void generateR();
	void generateSG(unsigned int parallelClass, ParallelClassIds& tmpSG, InputSet blockIntersection);
	void printResolution();
	void printM();
	vector< ParallelClassIds > generateResolutionIds();
	void printSG();
	void printResolutionIds();
	vector< NodeSet > generateNodeSGContain(int rank);
};


#endif
