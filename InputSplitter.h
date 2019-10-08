#ifndef _MR_INPUTSPLITTER
#define _MR_INPUTSPLITTER

#include <fstream>
#include "CodedConfiguration.h"

using namespace std;

class InputSplitter {
private:
	const CodedConfiguration* conf;
	unsigned long long int bufferSize;

public:
	InputSplitter();
	~InputSplitter();
	void setConfiguration( const CodedConfiguration* configuration ) { conf = configuration; }
	void splitInputFile();

private:
	void createSplit( ifstream &inputFile, int splitNumber, unsigned long long int startIndex, unsigned long long int endIndex, char* buff );
};



#endif
