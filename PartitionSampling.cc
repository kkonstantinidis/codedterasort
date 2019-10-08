#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <assert.h>
#include <string.h>
 
#include "PartitionSampling.h"
#include "Common.h"

using namespace std;

PartitionSampling::PartitionSampling()
{
	conf = NULL;
}

PartitionSampling::~PartitionSampling()
{

}

PartitionList* PartitionSampling::createPartitions()
{  
	if ( !conf ) {
		cout << "Configuration has not been provided.\n";
		assert( false );
	}

	ifstream inputFile( conf->getInputPath(), ios::in | ios::binary | ios::ate );
	if ( !inputFile.is_open() ) {
		cout << "Cannot open input file " << conf->getInputPath() << endl;
		assert( false );
	}

	PartitionList keyList;
	long unsigned int keySize = conf->getKeySize();
	unsigned long fileSize = inputFile.tellg();
	long unsigned int numSamples = min( conf->getNumSamples(), fileSize / conf->getLineSize() );


	// Sample keys
	inputFile.seekg( 0 );
	long unsigned int gapSkip = ( ( fileSize / conf->getLineSize() ) / numSamples - 1 ) * conf->getLineSize();
	for ( long unsigned int i = 0; i < numSamples; i++ ) {
		unsigned char *keyBuff = new unsigned char [ keySize + 1 ];
		if ( keyBuff == NULL ) {
			cout << "Cannot allocate memory to sample keys.\n";
			assert( false );
		}
		inputFile.read( (char *) keyBuff, keySize );            
		keyBuff[ keySize ] = '\0';
		keyList.push_back( keyBuff );
		inputFile.seekg( conf->getValueSize() + gapSkip , ios_base::cur );
	}
	inputFile.close();

	
	// Sort sampled keys
	stable_sort( keyList.begin(), keyList.end(), cmpKey );

	// Partition keys
	PartitionList* partitions = new PartitionList;
	long unsigned int numPartitions = conf->getNumReducer();
	long unsigned int sizePartition = round( numSamples / numPartitions );
	for ( unsigned long int i = 1; i < numPartitions; i++ ) {
		unsigned char *keyBuff = new unsigned char [ keySize + 1 ];
		if ( keyBuff == NULL ) {
			cout << "Cannot allocate memory for partition keys.\n";
			assert( false );
		}
		memcpy( keyBuff, keyList.at( i * sizePartition ), keySize + 1 );
		partitions->push_back( keyBuff );
	}

	// Clean up
	for ( auto k = keyList.begin(); k != keyList.end(); ++k ) {
		delete [] (*k);
	}

	return partitions;
}


bool PartitionSampling::cmpKey( const unsigned char* keyl, const unsigned char* keyr )
{
	for ( unsigned long i = 0; keyl[i] != '\0'; i++ ) {
		if ( keyl[ i ] < keyr[ i ] ) {
			return true;
		}
		else if ( keyl[ i ] > keyr[ i ] ) {
			return false;
		}
	}
	return true;
}
