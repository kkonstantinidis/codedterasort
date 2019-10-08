#ifndef _MR_PARTITIONSAMPLING
#define _MR_PARTITIONSAMPLING

#include <vector>
#include "CodedConfiguration.h"
#include "Common.h"

using namespace std;

class PartitionSampling {
private:
	const CodedConfiguration* conf;  

public:
	PartitionSampling();
	~PartitionSampling();
	void setConfiguration( const CodedConfiguration* configuration ) { conf = configuration; }
	PartitionList* createPartitions();

private:
	static bool cmpKey( const unsigned char* keyl, const unsigned char* keyr );
};

#endif
