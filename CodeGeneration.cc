#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>

#include "CodeGeneration.h"

//Needed for pow()
#include <cmath>

using namespace std;

CodeGeneration::CodeGeneration( int _N, int _K, int _k, int _q ): N( _N ), K( _K ), k( _k ), q( _q )
{
	generateR();
	
	//Test
	// printM();
	
	resolutionIds = generateResolutionIds();
	
	//Test
	// printResolutionIds();
	
	//Note the space allocation for the vector below
	vector< BlockId > tmpSG(k);
	
	InputSet blockIntersection;
	generateSG(0, tmpSG, blockIntersection);
	
	//Test 
	// printSG();
	
	//Assign an ID to each shuffling group, starting from 0
	for (unsigned int i = 0; i < shufflingGroups.size(); i++){
		SGIdMap[shufflingGroups[i]] = i; 
	}
	
	//Assign each slave a vector of the shuffling groups it participates into
	for (int rank = 1; rank <= K; rank++){
		NodeSGMap[rank] = generateNodeSGContain(rank);
	}
	
	//Test
	// printf("Some slave finished CodeGen\n");
	
}


SGId CodeGeneration::getSGId(NodeSet ns){
	auto it = SGIdMap.find(ns);
	if(it != SGIdMap.end()){
		return it->second;
	}else{
		cout << "Cannot find SubsetId\n";
		assert(false);
	}
}


vector< NodeSet > CodeGeneration::generateNodeSGContain(int rank){
	
	//The variable to be returned
	vector< NodeSet > ret;
	
	//Test (prints only the IDs)
	// printf("Slave %d is assigned the shuffling groups with IDs:{ ", rank);
	
	//Go through each shuffling group
	for (auto SGit = shufflingGroups.begin(); SGit != shufflingGroups.end(); ++SGit){
		
		//If you find current slave's rank, store the corresponding shuffling group.
		if ((*SGit).find(rank) != (*SGit).end()){
			
			//Test
			// printf("%d ", SGIdMap[*SGit]);
			
			ret.push_back(*SGit);
		}
		
	}
	
	//Test
	// printf("}\n");
	
	return ret;
	
}


//Test
void CodeGeneration::printResolutionIds(){
	for (int i = 0; i < k; i++){
		printf("Block IDs for parallel class %d: {", i);		
		for(auto j = resolutionIds.at(i).begin(); j != resolutionIds.at(i).end(); ++j){
			printf("%d ", *j);
		}
		printf("}\n");
	}
}


//Test
void CodeGeneration::printSG(){
	int ctr = 0;
	for (auto sg = shufflingGroups.begin(); sg != shufflingGroups.end(); ++sg){
		NodeSet ns = *sg;
		printf("Shuffling group %d contains the slaves: {", ctr);
		for(auto slave = ns.begin(); slave != ns.end(); ++slave){
			printf("%d ", *slave);
		}
		printf("}\n");
		ctr++;
	}
}


/*
It will generate all possible sets of k blocks from distinct parallel classes that satisfy Shuffling phase stage 1 (vector of sets of int). The integers will be the IDs of the blocks (ranks of 
the slaves). 
*/ 
void CodeGeneration::generateSG(unsigned int parallelClass, ParallelClassIds& tmpSG, InputSet blockIntersection){
	
	for (auto block = resolutionIds[parallelClass].begin(); block != resolutionIds[parallelClass].end(); ++block){

		//Store the current intersection in order to be able to restore it at the next iteration when you will have to decide on another block from the last parallel class
		InputSet interTmpO = blockIntersection;
		
		tmpSG[parallelClass] = *block;
		
		//Test
		// printf("Parallel class %d and current block: %d\n\n", parallelClass, *block);
		
		if (blockIntersection.size() == 0){
			
			//Test
			// printf("Empty intersection, inserting block %d\n", *block);
			
			//If the intersection is empty, initialize it
			blockIntersection = M[*block];
			
		}else{
			
			//Test
			// printf("Non-empty intersection, computing intersection with block %d\n", *block);
			
			//Else, compute its intersection with the block you just selected. The temporary variable is needed since you can't store the intersection directly to itself
			InputSet interTmpI;
			set_intersection(blockIntersection.begin(), blockIntersection.end(), M[*block].begin(), M[*block].end(), inserter(interTmpI, interTmpI.begin()));
			blockIntersection = interTmpI;
			
		}
		
		//Test
		// printf("The intersection contains the subfiles: {");
		// for (auto it = blockIntersection.begin(); it != blockIntersection.end() ; ++it){
			// printf("%d ", *it);
		// }
		// printf("}\n");
		
		if (parallelClass == resolutionIds.size() - 2) {
			
			//Now it's time to get the block from last parallel class and store the group. There are q-1 ways to do that.
			
			//Test
			// printf("Chosen blocks from k-1 first parallel classes and now from class %d\n", parallelClass+1);
			
			bool commonFileFound = false;
			for (auto lastBlock = resolutionIds[parallelClass+1].begin(); lastBlock !=  resolutionIds[parallelClass+1].end(); ++lastBlock){
				
				//This is done to avoid some intersection computations since the common subfile is only in one block, and as soon as we find we don't need to do other intersections
				InputSet lastBlockIntersection;
				if(!commonFileFound){ 
					set_intersection(blockIntersection.begin(), blockIntersection.end(), M[*lastBlock].begin(), M[*lastBlock].end(), 
					inserter(lastBlockIntersection, lastBlockIntersection.begin()));
				}
				
				//If the common subfile is here, do not add this block
				if (lastBlockIntersection.size() == 1){
					
					//Test
					// printf("Block %d is rejected\n", *lastBlock);
					
					commonFileFound = true;
					continue;
				}
				
				//Test
				// printf("Block %d has been chosen from last parallel class\n", *lastBlock);
				
				//This block can be added as the last one
				tmpSG[parallelClass+1] = *lastBlock;
				
				//Convert vector to set. Alternatively copy(tmpSG.begin(), tmpSG.end(), inserter(tmpSGS, tmpSGS.end()))
				NodeSet tmpSGS(tmpSG.begin(), tmpSG.end());
				
				//Test
				// printf("Shuffling group contains the slaves: {");
				// for (auto slave = tmpSGS.begin(); slave != tmpSGS.end(); ++slave){
					// printf("%d ", *slave);
				// }
				// printf("}\n");
				
				//Before storing the group...
				
				//...create all subgroups of k-1 users and map the supergroup to them
				int ctr = 0;
				for (unsigned int i = 0; i < tmpSGS.size(); i++){
					
					//Restore the subgroup to be the whole supergroup
					NodeSet tmpSGa(tmpSGS.begin(), tmpSGS.end());
					
					//Delete i-th slave from supergroup.
					auto tmpSGaIt = tmpSGa.begin();
					advance(tmpSGaIt, i);
					
					//Test
					// printf("Removing slave %d from group...\n", *tmpSGaIt);
					
					//Delete current element from group
					tmpSGa.erase(tmpSGaIt);

					//Find common subfile and map subgroup to their common subfile. We don't need to worry since the set is guaranteed to have only 1 element.
					auto tmpSGaIt2 = tmpSGa.begin();
					InputSet subgroupIntersection = InputSet(M[*tmpSGaIt2]);
					
					//Test
					// printf("The files of the 1st slave of the intersection are: {");
					// for (auto it = subgroupIntersection.begin(); it != subgroupIntersection.end() ; ++it){
						// printf("%d ", *it);
					// }
					// printf("}\n");
					
					tmpSGaIt2++;
					
					for (; tmpSGaIt2 != tmpSGa.end(); ++tmpSGaIt2){
						
						//Test
						// printf("Currently computing intersection with the files of slave %d...\n", *tmpSGaIt2);
						
						InputSet subgroupIntersectionTmp;
						set_intersection(subgroupIntersection.begin(), subgroupIntersection.end(), M[*tmpSGaIt2].begin(), M[*tmpSGaIt2].end(), inserter(subgroupIntersectionTmp, subgroupIntersectionTmp.begin()));

						//Important: You need to clear it first
						subgroupIntersection.clear();
						
						//Now copy...
						copy(subgroupIntersectionTmp.begin(), subgroupIntersectionTmp.end(), inserter(subgroupIntersection, subgroupIntersection.end()));
					}
					auto commonSubfileIt = subgroupIntersection.begin();
					commonSubfile[tmpSGa] = *commonSubfileIt;
					
					//Test
					// printf("The shuffling SUB-group contains the slaves: {");
					// for (auto slave = tmpSGa.begin(); slave != tmpSGa.end(); ++slave){
						// printf("%d ", *slave);
					// }
					// printf("} and their common subfile is %d\n", commonSubfile[tmpSGa]);
					
					ctr++; 
				}
				
				//The shuffling group set is ready
				shufflingGroups.push_back(tmpSGS);
				
			}
			
			
		}else{
			
			//Test
			// printf("Not done yet, proceed to class %d\n", parallelClass+1);
			
			//Select block from next parallel class
			generateSG(parallelClass+1, tmpSG, blockIntersection);
			
		}
		
		//Restore initial value for next iteration
		blockIntersection = interTmpO;
	}
	
}


vector< ParallelClassIds > CodeGeneration::generateResolutionIds(){
	
	//The variable to be returned
	vector< ParallelClassIds > ret(k, vector< BlockId >(q));
	
	for (int i = 0; i < k; i++){	
		for(int j = 0; j < q; j++){
			// printf("Current parallel class size %ld\n", ret.at(i).size());
			ret[i][j] = i*q+j+1;
		}
	}
	
	return ret;
}


//Test
void CodeGeneration::printM(){
	for (int i = 0; i < K; i++){
		printf("Slave %d has been assigned the subfiles: {", i+1);
		for (auto fileIterator = M[i+1].begin(); fileIterator != M[i+1].end(); fileIterator++){
			int fileId = *fileIterator;
			printf("%d ", fileId);
		}
		printf("}\n");
	}
}


void CodeGeneration::generateR(){

	//The variable to be returned
	vector< ParallelClass > ret;

	//The blocks will be stored here
	Block B[k][q];
	
	typedef vector<int> codeword;
	
	vector< codeword > codewords;

	//Generate 1st codeword
	codeword first_word;
	int slaveId = 1;
	for (int i = 0; i <= k-1; i++){
		first_word.push_back(0);

		//Add codeword to block
		B[i][0].push_back(0);
		
		//Assign file to slave
		M[slaveId].insert(1);
		
		//Assign slave to file
		FileNodes[1].insert(slaveId);
		
		//Test
		// printf("Slave %d has been assigned the file %d\n", slaveId, 1);
		
		slaveId += q;
	}
	codewords.push_back(first_word);

	int digits_sum = 0;

	//Generate the rest of the codewords
	for (int i = 1; i < (int)pow(q, k-1); i++){

		codeword current_word = codewords.at(i-1);

		int pos = k-2;
		bool done = false;
		while (pos >= 0){
			
			if (done == false){
				if (current_word.at(pos) < q-1){

					current_word.at(pos)++;

					digits_sum++;
					done = true;

				}else{

					digits_sum -= current_word.at(pos);
					current_word.at(pos) = 0;

				} 
			}

			//Add codeword to block
			B[pos][current_word.at(pos)].push_back(i);

			//Assign file to slave
			//For the row-wise mapping the resolution which is 2D to 1D 
			M[pos*q+current_word.at(pos)+1].insert(i+1);
			
			//Assign slave to file
			FileNodes[i+1].insert(pos*q+current_word.at(pos)+1);
			
			//Test
			// printf("Slave %d has been assigned the file %d\n", pos*q+current_word.at(pos)+1, i+1);
			
			pos--;

		}
		current_word.at(k-1) = digits_sum % q;

		//Add codeword to block
		B[k-1][current_word.at(k-1)].push_back(i);
		
		//Assign file to slave
		//For the row-wise mapping the resolution which is 2D to 1D 
		M[(k-1)*q+current_word.at(k-1)+1].insert(i+1);
		
		//Assign file to slave
		FileNodes[i+1].insert((k-1)*q+current_word.at(k-1)+1);

		//Test
		// printf("Slave %d has been assigned the file %d\n", (k-1)*q+current_word.at(k-1)+1, i+1);
		
		codewords.push_back(current_word); 

	}
	
	//Test
	// printf("\nThe blocks are: \n");
	// for (int i = 0; i <= k-1; i++){
		// for (int j = 0; j <= q-1; j++){
			// printf("Block[%d][%d]: { ", i+1, j);
			
			// for (unsigned int m = 0; m < B[i][j].size(); m++){
				// printf("%d ", B[i][j].at(m));
			// }
			// printf("}\n");
		// }	
	// }
	// printf("\n");
	
	//Test
	// int m = 1;
	// printf("\nThe parallel classes are: ");
	// for (int i = 0; i < k; i++){
		// for (int j = 0; j < q; j++){
			// if (j % q == 0){
				// printf("\nP_%d:\n", m);
			// }
			// printf("Block[%d][%d]: { ", i+1, j);

			// for (unsigned int m = 0; m < B[i][j].size(); m++){
				// printf("%d ", B[i][j].at(m));
			// }

			// printf("}\n");
			// if (j % q == 0){
				// m++;
			// }
		// }
	// }
	
	//Test
	// for (int i = 1; i <= (int)pow(q, k-1); i++){
		// printf("The file %d has been assigned to the following slaves: ", i);
		// printNodeSet(FileNodes[i]);
		// printf("\n");
	// }
	
}


void CodeGeneration::printNodeSet( NodeSet ns )
{
	cout << '{';
	for( auto nit = ns.begin(); nit != ns.end(); ++nit ) {
		cout << ' ' << *nit;
		if ( nit != --ns.end() ) {
			cout << ',';
		}
	}
	cout << " }";
}
