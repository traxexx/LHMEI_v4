#include <iostream>
#include <string>
#include "ComputeLHMEI.h"

void LHMEI_Help();
void LHMEI_Help_Detailed();

int main(int argc, char * argv[])
{
	if (argc <= 1) {
		LHMEI_Help();
		return 0;
	}
	
	std::string ArgString;
	std::string Dummies;
	std::string PATH = std::string("/net/wonderland/home/saichen/LHMEI_v4");
	
	ArgString += "-Win=600;-Step=100;-CtrlChr=20;-Simplify=1;-NonOffset=1;-MElist=" + PATH + "/refs/MobileElement.list;-MEcoord=";
	ArgString += PATH + "/refs/MobileElement.coord;-Mapper=bwa;-HetIndex=";
	ArgString += PATH + "/refs/hs37d5-chr20-MEI-slice.het-index;-Chr=-1;";
	
	Dummies = std::string("--verbose;--debug;--keepIntermediates;--includeSingleAnchor");
	
	std::string FirstArg = std::string(argv[1]);
	if (FirstArg.compare("Test") == 0) {  // test mode (intermediate file is default kept)
		std::cout << std::endl;
		std::cout << "  Running test mode..." << std::endl;
		std::cout << std::endl;
		argc--;
		argv++;
		ArgString += std::string("-Sample=TestSample;-Bam=") + PATH + "/usage_test/1612_test.bam;-WorkDir=";
		ArgString += PATH + "/usage_test/output;-ProgramDir=" + PATH;
		ArgString += ";-Mapper=/net/wonderland/home/mktrost/dev/gotcloud/bin/bwa-aln;";
	}
	else if (FirstArg.compare("-h") == 0) { // display detailed help info
		std::cout << std::endl;
		LHMEI_Help_Detailed();
		return 0;
	}
	else {  // normal mode
		ArgString += std::string("-Sample= ;-Bam= ;-WorkDir= ;-ProgramDir=");
		ArgString += PATH + (";-Mapper=/net/wonderland/home/mktrost/dev/gotcloud/bin/bwa-aln");
	}
	Options MainOptions( argc, argv, ArgString, Dummies );
	
	Options * ptrMainOptions = &MainOptions;
	ComputeLHMEI(ptrMainOptions);

	return 0;	
}


using std::cout;
using std::endl;

void LHMEI_Help ()
{
	cout << endl;
	cout << "LHMEI Options: [] is default" << endl;
	cout << "    Test:  run test commands. Override all options below except for --options." << endl;
	cout << endl;
	cout << "  Required:" << endl;
	cout << "    -Sample:  Sample name." << endl;
	cout << "    -Bam:  Raw bam files for LHMEI discovery. MUST be sorted by cooridinate." << endl;
	cout << "    -WorkDir:  Working directory (if not exist, LHMEI will create one)." << endl;
	cout << "    -ProgramDir: Directory of LHMEI." << endl;
	cout << endl;
}

void LHMEI_Help_Detailed()
{
	LHMEI_Help();
	cout << "  With defaults:" << endl;
	cout << "    -Mapper:  Mapper fore genearating -Bam. [bwa]" << endl;
	cout << "    -HetIndex:  Contrl stats of sliced MEI. [/LHMEI-dir/refs/hs37d5-chr20-MEI-slice.het-index]" << endl;
	cout << "    -MElist:  MEI Consensus sequence list. [refs/MobileElement.list]" << endl;
	cout << "    -MEcoord:  Genomic regions of MEs. To speed up. Recommended. [ /refs/MobileElement.coord ]" << endl;
	cout << "    --CtrlChr: Chr to remap [20]" << endl;
	cout << endl;
	cout << "  Optional:" << endl;
	cout << "    -Simplify: Simplification level. The higher the more time-saving but less sensitivity. [1]" << endl;
	cout << "    -Chr: Only run discovery on this chr. If set -1, do whole genome. [-1]" << endl;
	cout << "    --verbose: Print out user-defined options. [Off]" << endl;
	cout << "    --debug: Print out detailed warning message. [Off]" << endl;
	cout << "    --keepIntermediates: Do not remove intermediate files. [Off]" << endl;
	cout << endl;
}
