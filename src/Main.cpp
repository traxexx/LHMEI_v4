#include <iostream>
#include <string>
#include "ComputeLHMEI.h"
#include "Wrappers.h"

using std::string;

int main(int argc, char * argv[])
{
	if (argc <= 1) {
		DisplayUsageInfo();
		return 0;
	}
	
	string Path = GetExePath(); // secured last is '/'
	if ( Path.length() <= 4 ) {
		std::cerr << "ERROR: LHMEI-Discovery is not in $ProgramDir/bin/" << std::endl;
		exit(1);
	}
	Path = Path.substr(0, Path.size() - 4); // remove bin/
	string RefPath = Path + "refs/";
	
	std::string ArgString;
	std::string Dummies;
	
	ArgString = "-Win=600;-Step=100;-CtrlChr=20;-Simplify=1;-NonOffset=1;-Chr=-1;";
	ArgString += "-MElist=" + RefPath + "MobileElement.list;-MEcoord=" + RefPath + "MobileElement.coord;-HetIndex=" + RefPath + "hs37d5-chr20-MEI-slice.het-index;";
	ArgString += "-SliceFA=" + RefPath + "slice-chr20-hs37d5.fa;";
	ArgString += "-Mapper=/net/wonderland/home/mktrost/dev/gotcloud/bin/bwa-mem;";
	ArgString += "-refPrefix=refStats;-ReadLen=-1;-InsSize=-1;-MeiType=-1;-Depth=-1;";
	
	Dummies = std::string("--verbose;--debug;--keepIntermediates;--includeSingleAnchor;--pseudoChr;--printNonVariant;--printRefStats;--noCtrlVcf;--disableDPfilter");
	
	std::string FirstArg = std::string(argv[1]);
	if (FirstArg.compare("Test") == 0) {  // test mode (intermediate file is default kept)
		std::cout << std::endl;
		std::cout << "  Running test mode..." << std::endl;
		std::cout << std::endl;
		argc--;
		argv++;
		ArgString += "-Sample=TestSample;-Bam=" + Path + "usage_test/1612_test.bam;-WorkDir=" + Path + "usage_test/output";
	}
	else if (FirstArg.compare("-h") == 0) { // display detailed help info
		std::cout << std::endl;
		DisplayDetailedUsageInfo();
		return 0;
	}
	else if (FirstArg.compare("-hello") == 0) { // display debug usage info
		DisplayDetailedUsageInfo();
		DisplayDebugUsageInfo();
		return 0;
	}
	else {  // normal mode
		ArgString += std::string("-Sample=.;-Bam= ;-WorkDir= ;");
	}
	
// do the work	
	Options MainOptions( argc, argv, ArgString, Dummies );
	Options * ptrMainOptions = &MainOptions;
	ComputeLHMEI(ptrMainOptions);

	return 0;	
}


