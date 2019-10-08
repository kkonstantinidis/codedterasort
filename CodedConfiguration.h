#ifndef _CMR_CONFIGURATION//aa
#define _CMR_CONFIGURATION 

#include <math.h>

using namespace std;

class CodedConfiguration{

private:

	unsigned int numReducer;
	unsigned int numInput;  
	const char *inputPath;
	const char *outputPath;
	const char *partitionPath;
	unsigned long numSamples;
	unsigned int load;
	
	//For some factorization K = k*q
	unsigned int q;
	unsigned int k;

public:
	
	CodedConfiguration(){  
		numReducer = 6;  // K
		
		//Test
		q = 2;
		k = 3;
		
		inputPath = "/home/ubuntu/Input/Input10000";
		
		outputPath = "/home/ubuntu/Output/Output10000";
		
		partitionPath = "/home/ubuntu/Partition/Partition10000";
		
		numSamples = 100;    
		
		//Factorize K as K = k*q
		// q = 2;
		// for(; q <= floor(sqrt(numReducer)); q++){
			// if(numReducer%q == 0) break;
		// }	
		// k = numReducer/q;
		
		numInput = (int)pow(q, k-1);//N = q^(k-1)  
		load = k;//==r
	}
	~CodedConfiguration(){}
	const static unsigned int KEY_SIZE = 10;
	const static unsigned int VALUE_SIZE = 90;  
	unsigned int getNumReducer()const{return numReducer;}
	unsigned int getNumInput()const{return numInput;}  
	const char *getInputPath()const{return inputPath;}
	const char *getOutputPath()const{return outputPath;}
	const char *getPartitionPath()const{return partitionPath;}
	unsigned int getKeySize()const{return KEY_SIZE;}
	unsigned int getValueSize()const{return VALUE_SIZE;}
	unsigned int getLineSize()const{return KEY_SIZE + VALUE_SIZE;}
	unsigned long getNumSamples()const{return numSamples;}  
	unsigned int get_q()const{return q;}
	unsigned int get_k()const{return k;}
	unsigned int getLoad()const{return load;}
};

#endif
