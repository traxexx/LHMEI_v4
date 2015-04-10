#include "OriginalStats.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <utility>
#include "Utilities.h"
#include "GLs.h"
#include "Globals.h"

using std::cout;
using std::cerr;
using std::endl;
using std::to_string;
using std::stoi;
using std::ifstream;
using std::stringstream;
using std::getline;

// static members
const int OriginalStats::ClipStart = 2;
const int OriginalStats::DiscStart = 10;
const int OriginalStats::RawCellSize = 18;


// constructor: only set mei index
OriginalStats::OriginalStats( int mei_type, string & sample_name ):
	mei_index( mei_type ),
	SampleName( sample_name ),
	current_add_start( 0 )
{
	rawStats.clear();
	chrNameHash.clear();
	chrIndexVec.clear();
	if ( mei_index == 0 )
		mei_name = string("Alu");
	else if ( mei_index == 1 )
		mei_name = string("L1");
	else if ( mei_index == 2 )
		mei_name = string("SVA");
	else {
		cerr << "ERROR: mei_index = " << mei_index << ", should be >=0 & <=2 !" << endl;
		exit(1);
	}
}

// destructor
OriginalStats::~OriginalStats() {}


// add
void OriginalStats::Add( string current_chr, string & proper_name, string & disc_name )
{
	ifstream properFile;
	string proper_zero = string(proper_name) + ".0";
	properFile.open( proper_zero.c_str() );
	CheckInputFileStatus( properFile, proper_name.c_str() );
// add .0 first
	int current_chr_index = chrNameHash.size();
	if ( chrNameHash.find( current_chr ) != chrNameHash.end() ) {
		cerr << "ERROR: chr" << current_chr << " already exists!" << endl;
		exit(1);
	}
	chrNameHash[ current_chr ] = current_chr_index;
  // new cell for add
	RawCell * add_raw = new RawCell;
	add_raw->chr_index = current_chr_index;
	add_raw->win_index = 0;
	add_raw->counts.resize( RawCellSize, 0);
// for adding the next 2
	current_add_start = rawStats.size();
	string line;
	while( getline( properFile,line ) ) {
		stringstream ss;
		ss << line;
		string str_cord_index;
		getline( ss, str_cord_index, '\t');
		string str_proper;
		getline( ss, str_proper, '\t');
		string str_short;
		getline( ss, str_short, '\t');
		add_raw->win_index = stoi( str_cord_index );
		add_raw->counts[0] = stoi( str_proper);
		add_raw->counts[1] = stoi( str_short);
		rawStats.push_back( *add_raw );
	}
	delete add_raw;
	properFile.close();

// the proper-chr.mei_index
	string clip_name = string(proper_name) + "." + to_string(mei_index + 1);
	appendRawStats( clip_name, ClipStart );
// then disc
	string new_disc_name = disc_name + string(".") + to_string(mei_index);
	appendRawStats( new_disc_name, DiscStart );
}	

// append, used in add clip & disc
void OriginalStats::appendRawStats( string & rec_name, int base )
{
	ifstream recFile;
	recFile.open( rec_name.c_str() );
	CheckInputFileStatus( recFile, rec_name.c_str() );
	vector< RawCell >::iterator rs = rawStats.begin();
	rs += current_add_start;
	bool reach_end = 0;
	string line;
	while( getline( recFile, line ) ) {
		stringstream ss;
		ss << line;
		string str_cord_index;
		getline( ss, str_cord_index, '\t');
		int cord_index = stoi( str_cord_index ); // index of current window to add
		
		if ( cord_index > rs->win_index ) { // find next
			while( cord_index > rs->win_index ) {
				rs++;
				if (rs == rawStats.end()) {
					reach_end = 1;
					break;
				}
			}
			if (reach_end) // no more to add (no available proper any more )
				break;
		}
		// now cord_index <= rs->win_index
		
		if ( rs->win_index == cord_index) { // add
			for( int i = 0; i < 8; i++ ) {
				string field;
				getline( ss, field, '\t' );
				rs->counts[ i + base ] = stoi( field );
			}
			rs++;
			if ( rs == rawStats.end() )
				break; // all added
		}
		else { // cord < ref, no clip in current cord
			continue;
		}
	}
	recFile.close();
}


// construct genome link
void OriginalStats::ReOrganize()
{
// make chrIndexHash
	buildChrIndexVec();
	
// sort rawStats
	std::sort( rawStats.begin(), rawStats.end(), sortRawStats );

	vector<int> dupVec; // record dup
	setDupVecForRawStats( dupVec );
	int empty = 0; // first several empty cells
	vector<RawCell>::iterator rs = rawStats.begin();
	int sum = 0;
	for( vector<int>::iterator it = rs->counts.begin(); it != rs->counts.end(); it++ )
		sum += *it;
	if (sum == 0)
		empty += dupVec[0];
// pre-allocate merge data
	if (empty == 0)
		MergeData.resize( dupVec.size() );
	else
		MergeData.resize( dupVec.size() - 1 );
		
// raw->merge & set genome location map
	vector<RawCell>::iterator raw_ptr = rawStats.begin();
	vector<int>::iterator dup_ptr = dupVec.begin();
	MergeCellPtr merge_ptr = MergeData.begin();
	if (empty != 0) {
		raw_ptr += dupVec[0];
		dup_ptr++;
	}
  // start: clear rawStats when finish copy (to save memory)
  	GenomeLocationCell new_cell;
  	for( ; dup_ptr != dupVec.end(); dup_ptr++, merge_ptr++ ) {
  		merge_ptr->counts = std::move( raw_ptr->counts );
  		merge_ptr->dups = (*dup_ptr);
  		for( int ct = 0; ct < (*dup_ptr); ct++ ) {
  			new_cell.win_index = raw_ptr->win_index;
  			new_cell.ptr = merge_ptr;
  			GenomeLocationMap[ convertChrIndexToName( raw_ptr->chr_index ) ].push_back( new_cell );
  			raw_ptr++;
  		}
  	}

// sort genome location map by chr
	for( map< string, vector< GenomeLocationCell > >::iterator map_it = GenomeLocationMap.begin(); map_it != GenomeLocationMap.end(); map_it++ ) {
		std::sort( map_it->second.begin(), map_it->second.end(), sortGenomeLocations );
	}
// clear
	rawStats.clear();
}

// clear the cells with level < LEVEL
void OriginalStats::ClearUnderLevelMergeCells()
{
	int erase_cell = 0;
	for( MergeCellPtr mptr = MergeData.begin(); mptr != MergeData.end(); mptr++ ) {
		int n_support = getSumSupportDiscs( mptr->counts ) + getSumSupportUnmaps(mptr->counts) + getSumSupportClips(mptr->counts);
		if ( n_support < LEVEL ) {
			mptr->counts.clear();
			erase_cell++;
		}
	}
	if (DEBUG_MODE)
		cout << "Removed under level cell = " << erase_cell << endl;
}

// print out
void OriginalStats::PrintGLasVcf( string & vcf_name )
{
// do print. Also merge consecutive windows. Print out the one with highest GL
	ofstream outVcf;
	outVcf.open( vcf_name.c_str() );
	CheckOutFileStatus( outVcf, vcf_name.c_str() );
	outVcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << SampleName << endl;
	int Times = WIN / STEP; // how many over lap windows
	for( map< string, vector< GenomeLocationCell > >::iterator map_it = GenomeLocationMap.begin(); map_it != GenomeLocationMap.end(); map_it++ ) {
		// single chr
		string chr_name = map_it->first;
		GlcPtr anchor; // mark the anchor
		SingleCellPrint* infoPtr = NULL; // if NULL, anchor has not assigned yet
		for( GlcPtr cell_it = map_it->second.begin(); cell_it != map_it->second.end(); cell_it++ ) {
			if( cell_it->ptr->GL.size() != 3 )
				continue;
//if(infoPtr) printSingleMergeCell( outVcf, cell_it, chr_name, infoPtr );		
			if ( cell_it->ptr->GL[0] >= cell_it->ptr->GL[1] && cell_it->ptr->GL[0] >= cell_it->ptr->GL[2] ) // skip the window with no variant
				continue;
			
			if ( infoPtr == NULL ) { // new anchor
				anchor = cell_it;
				infoPtr = new SingleCellPrint( anchor->win_index, anchor->ptr->GL, anchor->ptr->counts );
			}
			else { // compare with previous anchor
				if ( cell_it->win_index - infoPtr->anchor_end > Times ) { // not consecutive. start new anchor. print old anchor
					printSingleMergeCell( outVcf, anchor, chr_name, infoPtr );
					delete infoPtr;
					anchor = cell_it;
					infoPtr = new SingleCellPrint( anchor->win_index, anchor->ptr->GL, anchor->ptr->counts );
				}
				else { // consecutive. compare with anchor
					infoPtr->anchor_end = cell_it->win_index;
					infoPtr->wcount++;
					bool use_anchor;
					bool rank_success = setAnchorRank( use_anchor, infoPtr->gq_peak, anchor->ptr, cell_it->ptr );
					if ( rank_success ) { // comparable
						if ( !use_anchor ) { // adjust info ptr with new anchor
							anchor = cell_it;
							infoPtr->UpdateCoordWithNew( anchor->win_index );
						}
						infoPtr->UpdateBothEnd( cell_it->ptr->counts );
					}
					else { // equal window, widen CI, do not change anchor
						infoPtr->UpdateWithEqual( cell_it->win_index, cell_it->ptr->counts );
					}
				}
			}
		}
	// print out last anchor
		if ( infoPtr != NULL )
			printSingleMergeCell( outVcf, anchor, chr_name, infoPtr );
	}
	outVcf.close();
}

/*** print sub function ***/

// if set successfully, return 1;
// compare. If Anchor is bigger (use anchor), theRank = 1;
// rule: dosage --> posterior-variant --> %(disc + clip + unmap) --> less %proper --> depth --> can't set
bool OriginalStats::setAnchorRank( bool & theRank, int & gq_peak, MergeCellPtr & Anchor, MergeCellPtr & NewPtr )
{
// update gq_peak first
	int anchor_gq_sig_higher = 0; // for comparing posterior next
	int new_gq = GetVariantQuality( NewPtr->GL );
	if ( new_gq > gq_peak ) {
		if ( new_gq - gq_peak >= 10 )
			anchor_gq_sig_higher = -1;
		gq_peak = new_gq;
	}
	else if ( new_gq < gq_peak ) {
		if ( gq_peak - new_gq >= 10 )
			anchor_gq_sig_higher = 1;
	}

/* compare dosgae
	int anchor_dosage = GetAlleleDosage( Anchor->GL );
	int new_dosage = GetAlleleDosage( NewPtr->GL );
	if ( anchor_dosage > new_dosage ) {
		theRank = 1;
		return 1;
	}
	else if (anchor_dosage < new_dosage) {
		theRank = 0;
		return 1;
	}
*/
// posterior variant
	if ( anchor_gq_sig_higher > 0 ) {
		theRank = 1;
		return 1;
	}
	else if ( anchor_gq_sig_higher < 0 ) {
		theRank = 0;
		return 1;
	}
/*
	int anchor_posterior = GetVariantPosterior( Anchor->GL );
	int new_posterior = GetVariantPosterior( NewPtr->GL );
	if ( anchor_posterior > new_posterior )
		theRank = 1;
	else if ( anchor_posterior < new_posterior )
		theRank = 0;
	return 1;
*/

// calculate depth here but do not use as comparison
	int anchor_dp = GetVecSum( Anchor->counts );
	int new_dp = GetVecSum( NewPtr->counts );
		
// %support
	int anchor_support_frac = GetSupportReadFraction( Anchor->counts, anchor_dp );
	int new_support_frac = GetSupportReadFraction( NewPtr->counts, new_dp );
	if ( anchor_support_frac > new_support_frac ) {
		theRank = 1;
		return 1;
	}
	else if ( anchor_support_frac < new_support_frac ) {
		theRank = 0;
		return 1;
	}

// %proper
	int anchor_proper = GetProperReadFraction( Anchor->counts, anchor_dp );
	int new_proper = GetProperReadFraction( NewPtr->counts, new_dp );
	if ( anchor_proper < new_proper ) {
		theRank = 1;
		return 1;
	}
	else if ( anchor_proper > new_proper ) {
		theRank = 0;
		return 1;
	}	

// depth
	if ( anchor_dp > new_dp ) {
		theRank = 1;
		return 1;
	}
	else if ( new_dp > anchor_dp ) {
		theRank = 0;
		return 1;		
	}

	return 0;
}



// print a single record. sub-function
void OriginalStats::printSingleMergeCell( ofstream & outVcf, GlcPtr & Ptr, string & chr_name, SingleCellPrint * infoPtr )
{	
	MergeCellPtr Merge = Ptr->ptr;
	outVcf << chr_name << "\t" << infoPtr->central << "\t.\t.\t<INS:ME:" << mei_name << ">\t";
	outVcf << infoPtr->gq_peak << "\tPASS\tSVTYPE=INS;END=";
	outVcf << infoPtr->var_end << ";CIPOS=";
	outVcf << -infoPtr->ci << "," << infoPtr->ci;
//	outVcf << ";CIEND=" << -infoPtr->ci << "," << infoPtr->ci;
	int clip_count = getSumSupportClips( Merge->counts );
	int disc_count = getSumSupportDiscs( Merge->counts);
	int unmap_count = getSumSupportUnmaps( Merge->counts );
	outVcf << ";BOTH_END=" << infoPtr->both_end << ";CLIP=" << clip_count << ";DISC=" << disc_count << ";UNMAP=" << unmap_count << ";WCOUNT=" << infoPtr->wcount;
	outVcf << "\tGT:DP:GQ:PL\t";
	string genotype = GetGenotype( Merge->GL );
	vector<int> PL;
	PL.clear();
	SetPLsFromGL( PL, Merge->GL );
	int depth = GetVecSum( Merge->counts );  // simpy add reads together
	int gt_quality = GetGenotypeQuality( Merge->GL );
	outVcf << genotype << ":" << depth << ":" << gt_quality << ":" << PL[0] << "," << PL[1] << "," << PL[2];
	outVcf << endl;
}


/**** inner functions ****/

// build chrIndexVec from chrNameHash for later use
void OriginalStats::buildChrIndexVec()
{
	int vec_size = chrNameHash.size();
	chrIndexVec.resize( vec_size, string("") );
	for( map< string, int >::iterator map_it = chrNameHash.begin(); map_it != chrNameHash.end(); map_it++ ) {
		if ( map_it->second >= vec_size ) {
			cerr << "ERROR: chr index " << map_it->second << " larger than chrNameHash size!" << endl;
			exit(1);
		}
		chrIndexVec[ map_it->second ] = map_it->first;
	}
// sanity check: if any chrIndexVec element is uninitialized
	for( vector< string >::iterator it = chrIndexVec.begin(); it != chrIndexVec.end(); it++ ) {
		if ( it->size() == 0 ) {
			cerr << "ERROR: index " << ( it - chrIndexVec.begin() ) << " in chrIndexVec is not initialized!" << endl;
			exit(1);
		}
	}
}

// set dup vec by counting dups in sorted rawStats (inner)
void OriginalStats::setDupVecForRawStats( vector<int> & dupVec )
{
	dupVec.clear();
	vector<RawCell>::iterator raw = rawStats.begin();
	vector<RawCell>::iterator prev = raw;
	raw++;
	int dup = 1;
	for( ; raw != rawStats.end(); raw++, prev++ ) {
		if ( prev->counts == raw->counts ) // it is ok since int == method is defined
			dup++;
		else {
			dupVec.push_back( dup );
			dup = 1;
		}
	}
	dupVec.push_back( dup );
}



// sort rawStats: counts > chr > win
bool OriginalStats::sortRawStats( RawCell x, RawCell y )
{
	if ( x.counts.size() != y.counts.size() ) {
		cerr << "ERROR: sortRawStats x y do not have same counts!" << endl;
		exit(1);
	}
	vector<int>::iterator x_it = x.counts.begin();
	vector<int>::iterator y_it = y.counts.begin();
	for( ; x_it != x.counts.end(); x_it++, y_it++ ) {
		if ( *x_it < *y_it )
			return 1;
		else if ( *x_it > *y_it )
			return 0;
	}
// all the same then coord
	if ( x.chr_index < y.chr_index )
		return 1;
	else if ( x.chr_index > y.chr_index )
		return 0;
	else { // compare win index
		if ( x.win_index < y.win_index )
			return 1;
		else if ( x.win_index > y.win_index )
			return 0;
		else {
			cerr << "ERROR: sortRawStats: Same location record appear twice at: chr" << x.chr_index << ": " << x.win_index << endl;
			exit(1);
		}
	}
}


// sort GenomeLocation in each chr
bool OriginalStats::sortGenomeLocations( GenomeLocationCell x, GenomeLocationCell y )
{
	if ( x.win_index < y.win_index )
		return 1;
	else if ( x.win_index > y.win_index )
		return 0;
	else {
		cerr << "ERROR: sortGenomeLocation: Same location record appear twice at: " << x.win_index << endl;
		exit(1);
	}
}


/*** other utility functions ****/

// convert chr name to index
int OriginalStats::convertChrNameToIndex( string chr_name )
{
	if ( chrNameHash.find( chr_name ) == chrNameHash.end() ) {
		cerr << "ERROR: Can't find " << chr_name << " in chrNameHash!" << endl;
		exit(1);
	}
	int chr_index = chrNameHash[ chr_name ];
	return chr_index;
}

// convert chr index back to name
string OriginalStats::convertChrIndexToName( int chr_index )
{
	if ( chr_index >= int( chrIndexVec.size() ) ) {
		cerr << "ERROR: Can't find index " << chr_index << " in chrIndexVec!" << endl;
		exit(1);
	}
	string chr_name = chrIndexVec[ chr_index ];
	return chr_name;
}
















