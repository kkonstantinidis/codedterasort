#ifndef _CMR_WORKER
#define _CMR_WORKER

#include <mpi.h>
#include <unordered_map>
#include <pthread.h>

#include "CodedConfiguration.h"
#include "CodeGeneration.h"
#include "Common.h"
#include "Utility.h"
#include "Trie.h"

using namespace std;

class CodedWorker
{
public:
	typedef unordered_map< SGId, MPI::Intracomm > MulticastGroupMap;
	typedef unordered_map< unsigned int, LineList* > PartitionCollection; // key = destID
	typedef unordered_map< unsigned int, PartitionCollection > InputPartitionCollection; // key = inputID
	typedef struct _DataChunk {
		unsigned char* data;
		unsigned long long size; // number of lines
	} DataChunk;
	typedef map< VpairList, vector< DataChunk > > DataPartMap;
	typedef unordered_map< SGId, DataPartMap > NodeSetDataPartMap; // [Encode/Decode]PreData
	typedef map< Vpair, unsigned long long > VpairSizeMap;
	typedef struct _MetaData {
		VpairList vpList;
		VpairSizeMap vpSize;
		unsigned int partNumber; // { 1, 2, ... }
		unsigned long long size; // number of lines
	} MetaData;
	typedef struct _EnData {
		vector< MetaData > metaList;
		unsigned char* data;	  // encoded chunk
		unsigned long long size;      // in number of lines
		unsigned char* serialMeta;
		unsigned long long metaSize;  // in number of bytes
	} EnData;
	typedef unordered_map< SGId, EnData > NodeSetEnDataMap;  // SendData
	typedef unordered_map< SGId, vector< EnData > > NodeSetVecEnDataMap;  // RecvData
	
private:
	MPI::Intracomm workerComm;
	PartitionList partitionList;
	CodeGeneration* cg;
	MulticastGroupMap multicastGroupMap;
	TrieNode* trie;
	InputPartitionCollection inputPartitionCollection;
	NodeSetEnDataMap encodeDataSend;
	NodeSetVecEnDataMap encodeDataRecv;
	NodeSet localLoadSet;
	LineList localList;
	clock_t memTime;

public: // Because of thread
	const CodedConfiguration* conf;
	unsigned int rank;  
	NodeSetDataPartMap encodePreData;
	NodeSetDataPartMap decodePreData;    

public:
	CodedWorker( unsigned int _rank ): rank( _rank ) {}
	~CodedWorker();
	void setWorkerComm( MPI::Intracomm& comm ) { workerComm = comm; }
	void run();

private:
	void genMulticastGroup();
	TrieNode* buildTrie( PartitionList* partitionList, int lower, int upper, unsigned char* prefix, int prefixSize, int maxDepth );
	void execMap();
	void execReduce();
	void execEncoding();
	void execShuffle();
	void execDecoding();
	void sendEncodeData( EnData& endata, MPI::Intracomm& comm );
	void recvEncodeData( SGId nsid, unsigned int actId, MPI::Intracomm& comm );
	void outputLocalList();
};


#endif
