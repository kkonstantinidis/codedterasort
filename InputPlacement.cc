#include <mpi.h>
#include <iostream>
#include <fstream>

#include "CodedConfiguration.h"
#include "CodeGeneration.h"

#define BUFF_SIZE 1000000   // 1MB
#define MAX_FILE_PATH 128

using namespace std;

int main( int argc, char* argv[] )
{
	CodedConfiguration codedConf;
	CodedConfiguration* conf;
	CodeGeneration* cg = NULL;

	// Initialize configuration, which is known to all nodes
	conf = &codedConf;
	cg = new CodeGeneration( conf->getNumInput(), conf->getNumReducer(), conf->get_k(), conf->get_q() );


	// Initialize OpenMPI
	MPI::Init();
	int nodeRank = MPI::COMM_WORLD.Get_rank();
	unsigned int nodeTotal = MPI::COMM_WORLD.Get_size();

	if( nodeTotal != conf->getNumReducer() + 1 ) {
		if( nodeRank == 0 ) {
			cout << "The number of workers mismatches the number of processes.\n";
		}
		MPI::Finalize();
		return 0;    
	}

	// Master or Workers
	if( nodeRank == 0 ) {
		ifstream inputFile( conf->getInputPath(), ios::in | ios::binary | ios::ate );
		if ( !inputFile.is_open() ) {
			cout << "Cannot open input file " << conf->getInputPath() << endl;
			assert( false );
		}
		cout << "inputFile " << conf->getInputPath() << " is open\n";
		
		unsigned int numInput = conf->getNumInput();
		unsigned int LineSize = conf->getLineSize();
		unsigned long long int fileSize = inputFile.tellg();
		unsigned long long int numLine = fileSize / LineSize;

		// Assume that the number of lines is higher than the number of workers
		assert( numInput <= numLine );
		assert( numInput > 1 );        

		// Split input file equally except the last split that possible includes extra lines
		unsigned long long int splitSize = ( numLine / numInput ) * LineSize;

		char* buff = new char[ BUFF_SIZE ];
		if ( buff == NULL ) {
			cout << "Cannot allocate buffer for splitting\n";
			assert( false );
		}

		// Send each input
		for ( unsigned int i = 0; i < (unsigned int) numInput; i++ ) {
			unsigned int fid = i + 1;
			unsigned long long startIndex = i * splitSize;
			unsigned long long endIndex = ( i != numInput - 1 ) ? ( i + 1 ) * splitSize - 1 : fileSize - 1;
			unsigned long long totalSize = endIndex - startIndex + 1;
			unsigned long long currIndex = startIndex;      
			unsigned long long copySize;
			MPI::Intracomm mgComm;      

			cout << "Sending file " << fid << endl;      
			// Initial transmission
			// create multicast domain
			//NodeSet& ns = cg->getNodeSetFromFileID( fid );	
			int color = 1;
			mgComm = MPI::COMM_WORLD.Split( color, nodeRank );
			mgComm.Bcast( &totalSize, 1, MPI::UNSIGNED_LONG_LONG, 0 );

			// read split and send to workers
			inputFile.seekg( startIndex );      
			while( currIndex <= endIndex ) {
				copySize = min( (unsigned long long) BUFF_SIZE, endIndex - currIndex + 1 );
				inputFile.read( buff, copySize );
				currIndex += copySize;

				// send to workers
				mgComm.Bcast( buff, copySize, MPI::CHAR, 0 );
			}

			// Free multicast group
			mgComm.Free();
		}

		// cleanup
		inputFile.close();
		delete [] buff;  
	}
	else {
		// Worker side (receive input files)
		char filePath[ MAX_FILE_PATH ];
		unsigned long long totalSize;
		char* buff = new char[ BUFF_SIZE ];
		unsigned long long currSize;
		unsigned long long recvSize;


		//Code
		unsigned int numInput = conf->getNumInput();
		for( unsigned int i = 0; i < numInput; i++ ) {
			unsigned int fid = i + 1;

			// create multicast domain
			NodeSet& ns = cg->getNodeSetFromFileID( fid );
			
			//Test
			if (nodeRank == 1){	
				printf("I am slave %d and the split %d belongs to the following slaves: \n", nodeRank, fid); 
				CodeGeneration::printNodeSet(ns);
				printf("\n");
			}
			
			int color = ( ns.find( nodeRank ) != ns.end() ) ? 1 : 0;
			MPI::Intracomm mgComm = MPI::COMM_WORLD.Split( color, nodeRank );
			if( color == 0 ) {
				// No file to receive with respect to this subset
				continue;
			}
			else {
				mgComm.Bcast( &totalSize, 1, MPI::UNSIGNED_LONG_LONG, 0 );
				sprintf( filePath, "%s_%d", conf->getInputPath(), fid - 1 );
				ofstream inputFile( filePath, ios::out | ios::binary | ios::trunc );
				if ( !inputFile.is_open() ) {
					cout << "Cannot open input file " << filePath << endl;
					assert( false );
				}
				currSize = 0;
				while( currSize < totalSize ) {
					recvSize = min( (unsigned long long) BUFF_SIZE, totalSize - currSize );
					mgComm.Bcast( buff, recvSize, MPI::CHAR, 0 );
					inputFile.write( buff, recvSize );
					currSize += recvSize;
				}
				inputFile.close();
			}
			mgComm.Free();
		}
		
		delete [] buff;
	}
	MPI::Finalize();
	return 0;
}
