#include <iostream>
#include <mpi.h>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <assert.h>
#include <algorithm>
#include <ctime>
#include <string.h>
#include <cstdint>

#include "CodedWorker.h"
#include "CodedConfiguration.h"
#include "Common.h"
#include "Utility.h"
#include "CodeGeneration.h"

//For dbg
#include <unistd.h>

using namespace std;

CodedWorker::~CodedWorker()
{
	for ( auto it = partitionList.begin(); it != partitionList.end(); ++it ) {
		delete [] *it;
	}

	// Delete from inputPartitionCollection
	InputSet inputSet = cg->getM( rank );
	for ( auto init = inputSet.begin(); init != inputSet.end(); init++ ) {
		unsigned int inputId = *init;
		LineList* list = inputPartitionCollection[ inputId ][ rank - 1 ];
		for ( auto lit = list->begin(); lit != list->end(); lit++ ) {
			delete [] *lit;
		}
		delete list;
	}  

	// Delete from encodePreData
	for ( auto it = encodePreData.begin(); it != encodePreData.end(); it++ ) {
		DataPartMap dp = it->second;
		for ( auto it2 = dp.begin(); it2 != dp.end(); it2++ ) {
			vector< DataChunk >& vdc = it2->second;
			for ( auto dcit = vdc.begin(); dcit != vdc.end(); dcit++ ) {
				delete [] dcit->data;
			}
		}
	}

	// Delete from localList
	for ( auto it = localList.begin(); it != localList.end(); ++it ) {
		delete [] *it;
	}

	delete trie;    
	delete cg;
	delete conf;
}

void CodedWorker::run()
{
	
	// RECEIVE CONFIGURATION FROM MASTER
	conf = new CodedConfiguration;
	MPI::COMM_WORLD.Bcast( (void*) conf, sizeof( CodedConfiguration ), MPI::CHAR, 0 );


	// RECEIVE PARTITIONS FROM MASTER
	for ( unsigned int i = 1; i < conf->getNumReducer(); i++ ) {
		unsigned char* buff = new unsigned char[ conf->getKeySize() + 1 ];
		MPI::COMM_WORLD.Bcast( buff, conf->getKeySize() + 1, MPI::UNSIGNED_CHAR, 0 );
		partitionList.push_back( buff );
	}

	clock_t time;
	double rTime;  
	
	// GENERATE CODING SCHEME AND MULTICAST GROUPS
	time = clock();
	// if (rank==15){
	cg = new CodeGeneration( conf->getNumInput(), conf->getNumReducer(), conf->get_k(), conf->get_q() );
	// }
	genMulticastGroup();
	// printf("I am slave %d and finished genMulticastGroup\n", rank);
	time = clock() - time;
	rTime = double( time ) / CLOCKS_PER_SEC;
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );    


	// {
	// int i = 0;
	// char hostname[256];
	// gethostname(hostname, sizeof(hostname));
	// printf("%d: attach %d\n", rank, getpid());
	// fflush(stdout);
	// while (0 == i)
	// sleep(5);
	// }
	
	// EXECUTE MAP PHASE
	time = clock();
	memTime = 0;
	execMap();
	time = clock() - time - memTime;
	rTime = double( time ) / CLOCKS_PER_SEC;  
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );

	//send map memory time
	rTime = double( memTime ) / CLOCKS_PER_SEC;
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );

	
	// EXECUTE ENCODING PHASE
	time = clock();
	memTime = 0;
	execEncoding();
	time = clock() - time - memTime;
	rTime = double( time ) / CLOCKS_PER_SEC;
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );      

	//send encode memory time
	rTime = double( memTime ) / CLOCKS_PER_SEC;
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );


	// SHUFFLING PHASE
	execShuffle();


	// EXECUTE DECODING PHASE
	time = clock();
	memTime = 0;
	execDecoding();
	time = clock() - time - memTime;
	rTime = double( time ) / CLOCKS_PER_SEC;
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );     

	//send decode memory time
	rTime = double( memTime ) / CLOCKS_PER_SEC;
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );
	

	// REDUCE PHASE
	time = clock();
	execReduce();
	time = clock() - time;
	rTime = double( time ) / CLOCKS_PER_SEC;  
	MPI::COMM_WORLD.Gather( &rTime, 1, MPI::DOUBLE, NULL, 1, MPI::DOUBLE, 0 );

	outputLocalList();
	
	//Test
	// printf("I am rank %d and I have finished everything\n", rank-1);
}


void CodedWorker::genMulticastGroup(){
	map< NodeSet, SGId > sgidmap = cg->getSGIdMap();
	for (auto nsit = sgidmap.begin(); nsit != sgidmap.end(); nsit++){
		NodeSet ns = nsit->first;
		SGId nsid = nsit->second;
		
		MPI::Intracomm mgComm;
		if (ns.find(rank) != ns.end()){
			
			//If I am in the group
			mgComm = workerComm.Split(1, rank);
			
			//Test
			// printf("I am slave %d and belong to the communicator of group ID %d\n", rank, nsid);
			
		}else{
			
			//If I am not in the group, give MPI_UNDEFINED as the color. This speedsup the partitipation check at the Shuffle phase
			mgComm = workerComm.Split(MPI_UNDEFINED, rank);
			
		}
		multicastGroupMap[nsid] = mgComm;
	}
}


void CodedWorker::execMap()
{
	//Test
	// clock_t time_map = clock();
	
	// Get the set of subfiles IDs that correspond to your block and you need to process
	InputSet inputSet = cg->getM( rank );

	//Test
	// if (rank == 1) cout << "getM took me " << setw(10) << clock() - time_map << endl;
	// time_map = clock();
	
	// Build trie
	unsigned char prefix[conf->getKeySize()];
	trie = buildTrie(&partitionList, 0, partitionList.size(), prefix, 0, 2);

	//Test
	// if (rank == 1) cout << "buildTrie took me " << setw(10) << clock() - time_map << endl;
	// time_map = clock();
	
	// Read input files and partition data
	for (auto init = inputSet.begin(); init != inputSet.end(); init++){
		unsigned int inputId = *init;
		
		//Test
		// printf("I am slave with rank %d trying to open subfile %d\n", rank, inputId);
		
		// Read input
		char filePath[ MAX_FILE_PATH ];
		sprintf( filePath, "%s_%d", conf->getInputPath(), inputId - 1 );
		ifstream inputFile( filePath, ios::in | ios::binary | ios::ate );
		if ( !inputFile.is_open() ) {
			cout << rank << ": Cannot open input file " << filePath << endl;
			assert( false );
		}
		
		unsigned long fileSize = inputFile.tellg();
		unsigned long int lineSize = conf->getLineSize();
		unsigned long int numLine = fileSize/lineSize;
		inputFile.seekg( 0, ios::beg );
		PartitionCollection& pc = inputPartitionCollection[ inputId ];
		
		// Create lists of lines
		for ( unsigned int i = 0; i < conf->getNumReducer(); i++ ){
			memTime -= clock();
			pc[ i ] = new LineList;
			memTime += clock();
		}
		
		// Partition data in the input file
		for ( unsigned long i = 0; i < numLine; i++ ){
			memTime -= clock();
			unsigned char* buff = new unsigned char[ lineSize ];
			memTime += clock();
			inputFile.read( (char *)buff, lineSize );

			unsigned int wid = trie->findPartition( buff );
			
			//Debug
			/* 			if (rank == 15 && inputId == 40 && (wid < 0 || wid > 15)){
				printf("Rank %d at line %lu of split ID %d which has %lu lines. The current line is %s and belongs to partition %d\n", rank, i, inputId, numLine, buff, wid);
				// int waiting = 1;
				// while(waiting == 1) sleep(5);
			} */
			
			pc[wid]->push_back( buff );   
			
		}

		//Test
		// if (rank == 1) cout << "Mapping file " << inputId << " took me " << setw(10) << clock() - time_map << endl;
		// time_map = clock();
		
		//Remove unnecessarily lists (partitions associated with the other nodes having the file)
		NodeSet fsIndex = cg->getNodeSetFromFileID( inputId );
		for ( unsigned int i = 0; i < conf->getNumReducer(); i++ ) {
			if( i + 1 != rank && fsIndex.find(i+1) != fsIndex.end() ) {
				LineList* list = pc[i];
				for ( auto lit = list->begin(); lit != list->end(); lit++ ) {
					delete [] *lit;
				}
				delete list;	
			}
		}

		//Test
		// if (rank == 1) cout << "Clearing useless maps for file " << inputId << " took me " << setw(10) << clock() - time_map << endl;
		// time_map = clock();
		
		inputFile.close();
	}
}


void CodedWorker::execEncoding()
{
	vector< NodeSet > SGvec = cg->getNodeSGContain( rank );
	unsigned lineSize = conf->getLineSize();
	
	for ( auto nsit = SGvec.begin(); nsit != SGvec.end(); nsit++ ){
		SGId nsid = cg->getSGId( *nsit );
		unsigned long long maxSize = 0;
		
		//Construct chunks of input from data with index ns\{q}
		for ( auto qit = nsit->begin(); qit != nsit->end(); qit++ ) {
			if( (unsigned int) *qit == rank ){
				continue;
			}
			int destId = *qit;      
			NodeSet inputIdx( *nsit );
			inputIdx.erase( destId );
			
			int fid = cg->getCommonSubfileId( inputIdx );
			VpairList vplist;      
			vplist.push_back(Vpair(destId, fid));//You might need to take care of destId during the shuffling/decoding
			
			unsigned int partitionId = destId - 1;
			
			LineList* ll = inputPartitionCollection[ fid ][ partitionId ];
			auto lit = ll->begin();
			unsigned int numPart = conf->getLoad()-1;//Note the change i.e. that we will create k-1 chunks
			unsigned long long chunkSize = ll->size()/numPart;//Number of lines (not bytes)
			
			//First chunk to second last chunk
			for ( unsigned int ci = 0; ci < numPart - 1; ci++ ) {
				memTime -= clock();
				unsigned char* chunk = new unsigned char[ chunkSize * lineSize ];
				memTime += clock();
				for ( unsigned long long j = 0; j < chunkSize; j++ ) {
					memcpy( chunk + j*lineSize, *lit, lineSize );
					lit++;
				}
				DataChunk dc;
				dc.data = chunk;
				dc.size = chunkSize;
				
				encodePreData[ nsid ][ vplist ].push_back( dc );
			}

			//Last chunk
			unsigned long long lastChunkSize = ll->size() - chunkSize * ( numPart - 1 ); 
			memTime -= clock();			
			unsigned char* chunk = new unsigned char[ lastChunkSize * lineSize ];
			memTime += clock();
			for ( unsigned long long j = 0; j < lastChunkSize; j++ ) {
				memcpy( chunk + j*lineSize, *lit, lineSize );
				lit++;
			}
			DataChunk dc;
			dc.data = chunk;
			dc.size = lastChunkSize;
			encodePreData[ nsid ][ vplist ].push_back( dc );
			
			//Determine associated chunk of a worker (order in ns)
			unsigned int rankChunk = 0; // in [ 0, ... , r - 1 ]
			for ( auto it = inputIdx.begin(); it != inputIdx.end(); it++ ) {
				if( (unsigned int) *it == rank ) {
					break;
				}
				rankChunk++;
			}
			maxSize = max( maxSize, encodePreData[ nsid ][ vplist ][ rankChunk ].size );

			//Remove unused intermediate data from Map
			for ( auto lit = ll->begin(); lit != ll->end(); lit++ ) {
				delete [] *lit;
			}
			delete ll;
		}
		
		// Initialize encode data
		memTime -= clock();
		encodeDataSend[ nsid ].data = new unsigned char[ maxSize * lineSize ](); // Initial it with 0
		memTime += clock();
		encodeDataSend[ nsid ].size = maxSize;
		unsigned char* data = encodeDataSend[nsid].data;
		
		// Encode Data
		for (auto qit = nsit->begin(); qit != nsit->end(); qit++){
			if ((unsigned int) *qit == rank){
				continue;
			}
			int destId = *qit;      
			NodeSet inputIdx(*nsit);
			inputIdx.erase(destId);
			unsigned long fid = cg->getCommonSubfileId(inputIdx);
			VpairList vplist;
			vplist.push_back(Vpair(destId, fid));//You might need to take care of destId during the shuffling/decoding

			//Determine associated chunk of a worker ( order in ns )
			unsigned int rankChunk = 0;//in [0, ... , r-1]
			for (auto it = inputIdx.begin(); it != inputIdx.end(); it++){
				if ((unsigned int) *it == rank){
					break;
				}
				rankChunk++;
			}
			
			//Start encoding
			unsigned char* predata = encodePreData[nsid][vplist][rankChunk].data;
			unsigned long long size = encodePreData[nsid][vplist][rankChunk].size;
			unsigned long long maxiter = size*lineSize/sizeof(uint32_t);
			for (unsigned long long i = 0; i < maxiter; i++){
				((uint32_t*)data)[i] ^= ((uint32_t*)predata)[i];
			}

			//Fill metadata
			MetaData md;
			md.vpList = vplist;
			md.vpSize[vplist[0]] = size;//Assume Eta = 1;
			md.partNumber = rankChunk + 1;
			md.size = size; 
			encodeDataSend[nsid].metaList.push_back(md);
		}
		
		// Serialize Metadata
		EnData& endata = encodeDataSend[ nsid ];
		unsigned int ms = 0;
		ms += sizeof( unsigned int ); // metaList.size()
		for ( unsigned int m = 0; m < endata.metaList.size(); m++ ) {
			ms += sizeof( unsigned int ); // vpList.size() //should be some pointer space
			ms += sizeof( int ) * 2 * endata.metaList[ m ].vpList.size(); // vpList
			ms += sizeof( unsigned int ); // vpSize.size() //should be some pointer space
			ms += ( sizeof( int ) * 2 + sizeof( unsigned long long ) ) * endata.metaList[ m ].vpSize.size(); // vpSize
			ms += sizeof( unsigned int ); // partNumber
			ms += sizeof( unsigned long long ); // size
		}
		encodeDataSend[ nsid ].metaSize = ms;

		memTime -= clock();
		unsigned char* mbuff = new unsigned char[ ms ];
		memTime += clock();
		unsigned char* p = mbuff;
		unsigned int metaSize = endata.metaList.size();
		memcpy( p, &metaSize, sizeof( unsigned int ) );
		p += sizeof( unsigned int );
		// meta data List
		for ( unsigned int m = 0; m < metaSize; m++ ) {
			MetaData mdata = endata.metaList[ m ];
			unsigned int numVp = mdata.vpList.size();
			memcpy( p, &numVp, sizeof( unsigned int ) );
			p += sizeof( unsigned int );
			// vpair List
			for ( unsigned int v = 0; v < numVp; v++ ) {//this loop has always 1 iteration, due to the fact that every vpList gets only 1 element 
				memcpy( p, &( mdata.vpList[ v ].first ), sizeof( int ) );
				p += sizeof( int );
				memcpy( p, &( mdata.vpList[ v ].second ), sizeof( int ) );
				p += sizeof( int );
			}
			// vpair size Map
			unsigned int numVps = mdata.vpSize.size();
			memcpy( p, &numVps, sizeof( unsigned int ) );
			p += sizeof( unsigned int );
			for ( auto vpsit = mdata.vpSize.begin(); vpsit != mdata.vpSize.end(); vpsit++ ) {
				Vpair vp = vpsit->first;
				unsigned long long size = vpsit->second;
				memcpy( p, &( vp.first ), sizeof( int ) );
				p += sizeof( int );
				memcpy( p, &( vp.second ), sizeof( int ) );
				p += sizeof( int );
				memcpy( p, &size, sizeof( unsigned long long ) );
				p += sizeof( unsigned long long );
			}
			memcpy( p, &( mdata.partNumber ), sizeof( unsigned int ) );
			p += sizeof( unsigned int );
			memcpy( p, &( mdata.size ), sizeof( unsigned long long ) );
			p += sizeof( unsigned long long );
		}
		encodeDataSend[ nsid ].serialMeta = mbuff;
	}
}


void CodedWorker::execShuffle(){
	// NODE-BY-NODE
	clock_t time;

	for (unsigned int activeId = 1; activeId <= conf->getNumReducer(); activeId++){
		unsigned long long tolSize;
		clock_t txTime;
		workerComm.Barrier();
		if (rank == activeId){
			time = clock();
			txTime = 0;
			tolSize = 0;
		}
		vector< NodeSet >& vset = cg->getNodeSGContain(activeId);

		for (auto nsit = vset.begin(); nsit != vset.end(); nsit++){

			NodeSet ns = *nsit;
			SGId nsid = cg->getSGId(ns);

			MPI::Intracomm mcComm = multicastGroupMap[nsid];
			
			//Ignore subset if I am not in that group
			if (mcComm == MPI::COMM_NULL){
				continue;
			}

			if (rank == activeId){
				txTime -= clock();
				sendEncodeData(encodeDataSend[nsid], mcComm);
				txTime += clock();
				EnData& endata = encodeDataSend[nsid];
				tolSize += (endata.size*conf->getLineSize()) + endata.metaSize + (2*sizeof(unsigned long long));
			}else{
				//Convert activeId to rootId of a particular multicast group
				unsigned int rootId = 0;
				for(auto nid = ns.begin(); nid != ns.end(); nid++){
					if((unsigned int)(*nid) == activeId){
						break;
					}
					rootId++;	  
				}
				recvEncodeData(nsid, rootId, mcComm);
			}
		}

		//Active node should stop timer here
		workerComm.Barrier();        
		if (rank == activeId){
			time = clock() - time;
			double rTime = double(time)/CLOCKS_PER_SEC;
			double txRate = (tolSize*8*1e-6)/(double(txTime)/CLOCKS_PER_SEC);
			MPI::COMM_WORLD.Send(&rTime, 1, MPI::DOUBLE, 0, 0);
			MPI::COMM_WORLD.Send(&txRate, 1, MPI::DOUBLE, 0, 0);      
		}
	}
}


void CodedWorker::execDecoding(){
	for (auto nsit = encodeDataRecv.begin(); nsit != encodeDataRecv.end(); nsit++){
		SGId nsid = nsit->first;
		vector< EnData >& endataList = nsit->second;
		for (auto eit = endataList.begin(); eit != endataList.end(); eit++){
			EnData& endata = *eit;
			unsigned char* cdData = endata.data;
			unsigned long long cdSize = endata.size;
			vector< MetaData >& metaList = endata.metaList;

			unsigned int numDecode = 0;
			//Decode per VpairList
			MetaData dcMeta;
			for (auto mit = metaList.begin(); mit != metaList.end(); mit++){
				MetaData& meta = *mit;
				if (encodePreData[nsid].find(meta.vpList) == encodePreData[nsid].end()){
					dcMeta = meta;
					//No original data for decoding;
					continue;
				}
				unsigned char* oData = encodePreData[nsid][meta.vpList][meta.partNumber-1].data;
				unsigned long long oSize = encodePreData[nsid][meta.vpList][meta.partNumber-1].size;
				unsigned long long maxByte = min(oSize, cdSize)*conf->getLineSize();
				unsigned long long maxIter = maxByte/sizeof(uint32_t);
				for(unsigned long long i = 0; i < maxIter; i++){
					((uint32_t*) cdData)[i] ^= ((uint32_t*) oData)[i];
				}
				numDecode++;
			}

			//sanity check
			if (numDecode != metaList.size()-1){
				cout << rank << ": Decode error " << numDecode << '/' << metaList.size() - 1 << endl;
				assert(numDecode != metaList.size()-1);
			}

			if (decodePreData[nsid][dcMeta.vpList].empty()){
				for (unsigned int i = 0; i < conf->getLoad(); i++){
					decodePreData[nsid][dcMeta.vpList].push_back(DataChunk());
				}
			}

			decodePreData[nsid][dcMeta.vpList][dcMeta.partNumber-1].data = cdData;
			decodePreData[nsid][dcMeta.vpList][dcMeta.partNumber-1].size = dcMeta.size;
		}
	}

	unsigned int partitionId = rank-1;
	unsigned int lineSize = conf->getLineSize();  

	//Get partitioned data from input files, already stored in memory.
	InputSet inputSet = cg->getM(rank);
	for (auto init = inputSet.begin(); init != inputSet.end(); init++){
		unsigned int inputId = *init;
		LineList* ll = inputPartitionCollection[inputId][partitionId];
		//copy line by line
		for (auto lit = ll->begin(); lit != ll->end(); lit++) {
			memTime -= clock();
			unsigned char* buff = new unsigned char[lineSize];
			memTime += clock();
			memcpy(buff, *lit, lineSize);
			localList.push_back(buff);
		}
		localLoadSet.insert(inputId);
	}

	//Get partitioned data from other workers
	for (auto nvit = decodePreData.begin(); nvit != decodePreData.end(); nvit++){
		DataPartMap& dpMap = nvit->second;
		for (auto vvit = dpMap.begin(); vvit != dpMap.end(); vvit++){
			VpairList vplist = vvit->first;
			vector< DataChunk > vdc = vvit->second;
			//Add inputId to localLoadSet
			for (auto vpit = vplist.begin(); vpit != vplist.end(); vpit++){
				localLoadSet.insert(vpit->second);
			}
			//Add data from each part to locallist
			for (auto dcit = vdc.begin(); dcit != vdc.end(); dcit++){
				unsigned char* data = dcit->data;
				for (unsigned long long i = 0; i < dcit->size; i++){
					memTime -= clock();
					unsigned char* buff = new unsigned char[lineSize];
					memTime += clock();
					memcpy(buff, data + i*lineSize, lineSize);
					localList.push_back(buff);	  	  
				}
				delete [] dcit->data;
			}
		}
	}

	if (localLoadSet.size() != conf->getNumInput()){
		cout << rank << ": Only have paritioned data from ";
		CodeGeneration::printNodeSet(localLoadSet);
		cout << endl;
		assert(false);
	}    
}


void CodedWorker::execReduce(){  
	sort(localList.begin(), localList.end(), Sorter(conf->getKeySize()));
}


void CodedWorker::sendEncodeData( CodedWorker::EnData& endata, MPI::Intracomm& comm )
{
	// Send actual data
	unsigned lineSize = conf->getLineSize();  
	int rootId = comm.Get_rank();
	comm.Bcast( &( endata.size ), 1, MPI::UNSIGNED_LONG_LONG, rootId );
	comm.Bcast( endata.data, endata.size*lineSize, MPI::UNSIGNED_CHAR, rootId );
	delete [] endata.data;

	// Send serialized meta data
	comm.Bcast( &( endata.metaSize ), 1, MPI::UNSIGNED_LONG_LONG, rootId ); 
	comm.Bcast( endata.serialMeta, endata.metaSize, MPI::UNSIGNED_CHAR, rootId );
	delete [] endata.serialMeta;
}


void CodedWorker::recvEncodeData( SGId nsid, unsigned int rootId, MPI::Intracomm& comm )
{
	EnData endata;
	unsigned lineSize = conf->getLineSize();

	// Receive actual data
	comm.Bcast( &( endata.size ), 1, MPI::UNSIGNED_LONG_LONG, rootId );
	endata.data = new unsigned char[ endata.size * lineSize ];
	comm.Bcast( endata.data, endata.size*lineSize, MPI::UNSIGNED_CHAR, rootId );

	// Receive serialized meta data
	comm.Bcast( &( endata.metaSize ), 1, MPI::UNSIGNED_LONG_LONG, rootId );
	endata.serialMeta = new unsigned char[ endata.metaSize ];
	comm.Bcast( ( unsigned char* ) endata.serialMeta, endata.metaSize, MPI::UNSIGNED_CHAR, rootId );

	// De-serialized meta data
	unsigned char* p = endata.serialMeta;
	unsigned int metaNum;
	memcpy( &metaNum, p, sizeof( unsigned int ) );
	p += sizeof( unsigned int );
	// meta data List  
	for ( unsigned int m = 0; m < metaNum; m++ ) {
		MetaData mdata;
		// vpair List
		unsigned int numVp;
		memcpy( &numVp, p, sizeof( unsigned int ) );
		p += sizeof( unsigned int );
		for ( unsigned int v = 0; v < numVp; v++ ) {
			Vpair vp;
			memcpy( &( vp.first ), p, sizeof( int ) );
			p += sizeof( int );
			memcpy( &( vp.second ), p, sizeof( int ) );
			p += sizeof( int );
			mdata.vpList.push_back( vp );
		}
		// VpairSize Map
		unsigned int numVps;
		memcpy( &numVps, p, sizeof( unsigned int ) );
		p += sizeof( unsigned int );
		for ( unsigned int vs = 0; vs < numVps; vs++ ) {
			Vpair vp;
			unsigned long long size;
			memcpy( &( vp.first ), p, sizeof( int ) );
			p += sizeof( int );
			memcpy( &( vp.second ), p, sizeof( int ) );
			p += sizeof( int );
			memcpy( &size, p, sizeof( unsigned long long ) );
			p += sizeof( unsigned long long );
			mdata.vpSize[ vp ] = size;
		}
		memcpy( &( mdata.partNumber ), p, sizeof( unsigned int ) );
		p += sizeof( unsigned int );
		memcpy( &( mdata.size ), p, sizeof( unsigned long long ) );
		p += sizeof( unsigned long long );
		endata.metaList.push_back( mdata );
	}
	delete [] endata.serialMeta;

	//Serial decoder
	encodeDataRecv[ nsid ].push_back( endata );
}


void CodedWorker::outputLocalList()
{
	char buff[ MAX_FILE_PATH ];
	sprintf( buff, "%s_%u", conf->getOutputPath(), rank - 1 );
	ofstream outputFile( buff, ios::out | ios::binary | ios::trunc );
	for ( auto it = localList.begin(); it != localList.end(); ++it ) {
		outputFile.write( ( char* ) *it, conf->getLineSize() );
	}
	outputFile.close();
	//cout << rank << ": outputFile " << buff << " is saved.\n";
}


TrieNode* CodedWorker::buildTrie( PartitionList* partitionList, int lower, int upper, unsigned char* prefix, int prefixSize, int maxDepth )
{
	if ( prefixSize >= maxDepth || lower == upper ) {
		return new LeafTrieNode( prefixSize, partitionList, lower, upper );
	}
	InnerTrieNode* result = new InnerTrieNode( prefixSize );
	int curr = lower;
	for ( unsigned char ch = 0; ch < 255; ch++ ) {
		prefix[ prefixSize ] = ch;
		lower = curr;
		while( curr < upper ) {
			if( cmpKey( prefix, partitionList->at( curr ), prefixSize + 1 ) ) {
				break;
			}
			curr++;
		}
		result->setChild( ch, buildTrie( partitionList, lower, curr, prefix, prefixSize + 1, maxDepth ) );
	}
	prefix[ prefixSize ] = 255;
	result->setChild( 255, buildTrie( partitionList, curr, upper, prefix, prefixSize + 1, maxDepth ) );
	return result;
}
