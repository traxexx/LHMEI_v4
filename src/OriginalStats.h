#ifndef ORIGINALSTATS_H
#define ORIGINALSTATS_H

#include "AnchorRank.h" // merge window when print
#include "SamFile.h"
// processing stats that are directly read from stat file

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iterator>

using std::vector;
using std::map;
using std::string;
using std::ofstream;

typedef struct {
	int chr_index;
	int win_index;
	vector<int> counts;
} RawCell;


typedef struct {
	int dups;
	vector<int> counts;
	vector<float> GL;
} MergeCell;

typedef vector<MergeCell>::iterator MergeCellPtr;

typedef struct {
	MergeCellPtr ptr;
	int win_index;
} GenomeLocationCell;

typedef vector< GenomeLocationCell >::iterator GlcPtr;

class OriginalStats
{
  public:
  	OriginalStats( int mei_type, string & sample_name );
  	~OriginalStats();
  	
	bool Add( string current_chr, string & proper_name, string & disc_name );
	void ReOrganize();
	void ClearUnderLevelMergeCells();
	void PrintGLasVcf( string & vcf_name, string & bam_name, string & ref_fasta ); // print to vcf
// final data structure
	vector< MergeCell > MergeData;
	map< string, vector< GenomeLocationCell > > GenomeLocationMap;
	vector< MergeCell > SpecialMergeData; // for adjustment
	
  private:
  	void printSingleMergeCell( ofstream & outVcf, GlcPtr & Ptr, string & chr_name, SingleCellPrint * infoPtr );
  	bool setAnchorRank( bool & theRank, int & gq_peak, MergeCellPtr & Anchor, MergeCellPtr & NewPtr );
  	
  	void appendRawStats( string & rec_name, int base );
  	void buildChrIndexVec();
  	void setDupVecForRawStats( vector<int> & dupVec );
  	static bool sortRawStats( RawCell x, RawCell y );
  	static bool sortGenomeLocations( GenomeLocationCell x, GenomeLocationCell y );
  	int convertChrNameToIndex( string chr_name );
  	string convertChrIndexToName( int chr_index );
  	int getBreakPointAndCI( string & chr_name, int & center, int & event_end, int & ci_low, int & ci_high ); // get break point from bam
  
 // constants 
  	const int ClipStart; // vector start of clip
  	const int DiscStart;
  	const int RawCellSize;
  	const float MinDepth; // for depth filter when print vcf
  	const float MaxDepth;
  	const int mei_index; // current mei type
  	string SampleName;
  	int current_add_start; // current start for adding rawStats
	vector<RawCell> rawStats; // start from rawStats end

// bam file header
	SamFile currentSam;
	SamFileHeader currentSamHeader;

	
// for name conversion
	map< string, int > chrNameHash;
	vector< string > chrIndexVec;
	string mei_name;
};

#endif