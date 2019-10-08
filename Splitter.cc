#include <iostream>

#include "CodedConfiguration.h"
#include "InputSplitter.h"

using namespace std;

int main( int argc, char* argv[] )
{
	CodedConfiguration codedConf;
	CodedConfiguration* conf;

	conf = &codedConf;
	cout << "Split data specified in CODED configuration\n";

	// SPLIT INPUT FILE TO N (ROUGHLY) EQUALLY FILES, WHERE N IS THE NUMBER OF WORKERS.
	cout << ": ----------\n";
	cout << ": Split input file\n";  
	InputSplitter inputSplitter;

	inputSplitter.setConfiguration( conf );
	inputSplitter.splitInputFile();
	cout << ":Done\n\n";

	return 0;
}


