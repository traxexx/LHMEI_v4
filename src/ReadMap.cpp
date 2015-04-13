#include "ReadMap.h"
#include "Utilities.h" // CheckFileStatus
#include "bamUtility.h"
#include "MeiCoord.h"
#include "Globals.h" // debug mode
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <algorithm> // find
#include <fstream>
#include <sstream>
#include <iterator>

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::cout;
using std::endl;

ReadMap::ReadMap( int winlen, int steplen, int avr_read_len, int avr_ins_size, const char* ref_chr, const char* mei_list ):
	REF_CHR(ref_chr),
	MinQuality(10),
	MinClip(20),
	ShortClip(10),
	win_len(winlen),
	step_len(steplen),
	ReadMapCellSize(50),
	AvrReadLen(avr_read_len),
	MinReadLen( avr_read_len/2 ),
	AvrInsSize( avr_ins_size )
{
	initializeMeiSeqRef( mei_list );
}

ReadMap::~ReadMap() {}


// for re-map
void ReadMap::initializeMeiSeqRef( const char* mei_list )
{
// initialize mei ref seq from mei_list
	ifstream mei_list_file;
	mei_list_file.open(mei_list);
	CheckInputFileStatus( mei_list_file, mei_list );
	std::string alu_fasta;
	std::string line_fasta;
	std::string sva_fasta;
	std::string line;
	while(std::getline(mei_list_file, line)) {
		std::stringstream ss; ss << line;
		std::string field;
		std::getline(ss, field, '\t');
		std::string mei_field = field;
		if (mei_field.compare("ALU") == 0) {
			std::getline(ss, field, '\t');
			alu_fasta = field;
		}
		else if (mei_field.compare("L1") == 0) {
			std::getline(ss, field, '\t');
			line_fasta = field;		
		}
		else if (mei_field.compare("SVA") == 0) {
			std::getline(ss, field, '\t');
			sva_fasta = field;		
		}
		else {
			std::cerr << "ERROR: Unidentified MEI name: " << mei_field << std::endl; exit(2);
		}
	}
	mei_list_file.close();
	if (alu_fasta.length() == 0 || line_fasta.length() == 0 || sva_fasta.length() == 0) {
		std::cerr << "ERROR: Invalid MEI fasta. Check " << mei_list << std::endl; exit(2);
	}
	const float map_ratio = 0.8; const float mismatch_ratio = 0.1;
	REF_SEQ.resize(3);
	REF_SEQ[0] = new RefSeq;
	REF_SEQ[0]->ReadSeq(alu_fasta.c_str(), map_ratio, mismatch_ratio);
	REF_SEQ[0]->addPolyAtail();	
	REF_SEQ[1] = new RefSeq;
	REF_SEQ[1]->ReadSeq(line_fasta.c_str(), map_ratio, mismatch_ratio);
	REF_SEQ[1]->addPolyAtail();
	REF_SEQ[2] = new RefSeq;
	REF_SEQ[2]->ReadSeq(sva_fasta.c_str(), map_ratio, mismatch_ratio);
	REF_SEQ[2]->addPolyAtail();
}


// here focus_chr = REF_CHR
// LEVEL = -1   
void ReadMap::SetMapFromCtrlBam( const char* ctrl_bam, const char* output_prefix, const char* mei_coord_list )
{
// file check
	SamFile samIn;
	samIn.OpenForRead(ctrl_bam);
	if ( !samIn.IsOpen() ) {
		cerr << "ERROR: can't open ctrl bam " << ctrl_bam << endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool ctrl_header_status = samIn.ReadHeader(samHeader);
	if (!ctrl_header_status) {
		cerr << "ERROR: Fail to read header: " << ctrl_bam << endl;
		exit(1);
	}
// bai status
	string bai = string(ctrl_bam) + ".bai";
	bool bai_status = samIn.ReadBamIndex( bai.c_str() );
	if ( !bai_status ) {
		cerr << "ERROR: ctrl bam not indexed!" << endl;
		exit(1);
	}

// output discBam	
	string disc_name = string(output_prefix) + "ctrl-disc.sam";
	SamFile discSam;
	discSam.OpenForWrite( disc_name.c_str() );	
	if ( !discSam.IsOpen() ) {
		cerr << "ERROR: can't open disc ctrl bam " << disc_name << endl;
		exit(1);
	}
	discSam.WriteHeader( samHeader);

// print disc sam
	MeiCoord * mei_coord = new MeiCoord( mei_coord_list );
	QC * DiscQC = new QC( AvrReadLen, AvrInsSize );
	for( int rec = 0; rec < samHeader.getNumSQs(); rec++ ) {
		String val = samHeader.getReferenceLabel( rec );
		string current_chr = string(val.c_str());
		bool section_status = samIn.SetReadSection( current_chr.c_str() );
		if (!section_status) {
			std::cerr << "ERROR: Unable to set read section to chr " << current_chr << std::endl;
			exit(1);
		}
		mei_coord->ResetVecPtr( current_chr );
		SamRecord sam_rec;
		if ( current_chr.compare(REF_CHR) != 0 ) { // non-ref chr: only print mate mapped disc
			while(samIn.ReadRecord(samHeader, sam_rec)) {
				if ( sam_rec.getFlag() & 0x2 )
					continue;
				if ( sam_rec.getFlag() & 0x8 )
					continue;
				bool qc_pass = DiscQC->PassQC( sam_rec );
				if ( !qc_pass )
					continue;
				if ( sam_rec.getFlag() & 0x4 ) {
					string mate_name = sam_rec.getMateReferenceName();
					if ( mate_name.compare( REF_CHR ) == 0 ) {
						simplifySamRec(sam_rec);
						discSam.WriteRecord(samHeader, sam_rec);
					}
					continue;
				}
				bool disc_pass = DiscQC->DiscSamPass( sam_rec );
				if ( !disc_pass )
					continue;
				mei_coord->SetEMtag( sam_rec );
				discSam.WriteRecord(samHeader, sam_rec);
			}
		}
		else { // ref
			while(samIn.ReadRecord(samHeader, sam_rec)) { // ref chr: print mapped disc
				if ( sam_rec.getFlag() & 0x2 )
					continue;
				if ( sam_rec.getFlag() & 0x4 ) {
					simplifySamRec(sam_rec);
					discSam.WriteRecord(samHeader, sam_rec);
					continue;
				}
				bool qc_pass = DiscQC->PassQC( sam_rec );
				if ( !qc_pass )
					continue;
				if ( sam_rec.getFlag() & 0x8 ) {
					if ( sam_rec.getMapQuality() >= MinQuality ) {
						mei_coord->SetEMtag( sam_rec );
						simplifySamRec(sam_rec);
						discSam.WriteRecord(samHeader, sam_rec);
					}
					continue;
				}
				bool disc_pass = DiscQC->DiscSamPass( sam_rec );
				if ( !disc_pass )
					continue;
				mei_coord->SetEMtag( sam_rec );
				discSam.WriteRecord(samHeader, sam_rec);
			}
		}	
	}
	
	discSam.Close();
	delete DiscQC;
	delete mei_coord;

// set candidate
	discCandidates.clear();
	vector< std::pair<int, int> > candidate_region;
	candidate_region.resize(1);
	candidate_region[0].first = 30000;
	int ref_chr_length = atoi( samHeader.getSQTagValue("LN", REF_CHR) );
	candidate_region[0].second = ref_chr_length - 30000;
	discCandidates[REF_CHR] = candidate_region;

// read proper by section
	QC * ProperQC = new QC( AvrReadLen, AvrInsSize );
	
	string string_ref_chr_name = string(REF_CHR);
	properMap.clear();
	properMap.resize( (ref_chr_length - win_len) / step_len + 1 );
	processProperReadsBySection(samIn, samHeader, string_ref_chr_name, 0, candidate_region[0].second, ProperQC); // -30000 but it's ok
	ProperQC->PrintQCsummary();
	delete ProperQC;
	samIn.Close();
	string ctrl_proper_name = string(output_prefix) + "proper-" + REF_CHR;
	printProperMapBySection( ctrl_proper_name.c_str() );

// sort disc
	string nsort_disc_bam = SortBamByName( disc_name.c_str() );

// process disc
	if (DEBUG_MODE) {
		cout << "Processing disc reads in REF_CHR: " << REF_CHR << endl;
		cout << endl;
	}
	processDiscReads( nsort_disc_bam.c_str() );
	printDiscMap( output_prefix );
}


// set data from target bam
void ReadMap::SetMapFromBam 
	(const char* bam, const char * disc_bam, const char* output_prefix, const char * focus_chr, const char * level_file_prefix)
{
	cout << "Working on original bam..." << endl;
	SamFile samIn;
	samIn.OpenForRead(bam);
	if ( !samIn.IsOpen() ) {
		std::cerr << "ERROR: Unable to open input bam: " << bam << std::endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool bam_header_status = samIn.ReadHeader(samHeader);
	if (!bam_header_status) {
		cerr << "ERROR: Fail to read header: " << bam << endl;
		exit(1);
	}
	
	// check if exist .bai
	string bai_name = string( bam ) + ".bai";
	ifstream bai;
	bai.open(bai_name.c_str());
	if (!bai.good()) {
		cerr << "ERROR: " << bai_name << " do not exists. Please use samtools index [bam] to generate it." << endl;
		exit(1);
	}
	bai.close();
	
	// set bai
	bool bai_status = samIn.ReadBamIndex( bai_name.c_str() );
	if (!bai_status) {
		cerr << "ERROR: Invalid bam index: " << bai_name << endl;
		exit(1);
	}
	
	// check if focus_chr exist
	if ( strlen(focus_chr) > 0 ) {
		bool focus_exist = ExistChrInBam( samHeader, focus_chr );
		if (!focus_exist) {
			cerr << "ERROR: " << focus_chr << " do not exist in bam!" << endl;
			exit(1);
		}
	}
		
	// process proper reads in bam by section
	QC * ProperQC = new QC( AvrReadLen, AvrInsSize );
	for(int i=0; i< samHeader.getNumSQs(); i++) {  
		String chr_name_str = samHeader.getReferenceLabel(i);
		string chr_name = string(chr_name_str.c_str());
		if ( strlen(focus_chr) > 0 && chr_name.compare(focus_chr) != 0 )
			continue;
		vector< std::pair<int, int> > candidate_region;
		int chr_length = atoi( samHeader.getSQTagValue("LN", chr_name.c_str()) );
		string level_name = string(level_file_prefix) + "." + chr_name;
		setCandidateRegion( candidate_region, chr_length, level_name.c_str() );
	// here, add for activating disc map
		discCandidates[chr_name] = candidate_region;
		if (candidate_region.size() == 0) {
			std::cerr << "Warning: no available candidate regions in " << chr_name << " under level " << LEVEL << ", skip this chr!" << endl;
			continue;
		}
	// process by chr section
		properMap.clear();
		properMap.resize( chr_length / step_len + 1 );
		time_t raw_time;
		time(&raw_time);
		if (chr_name.size() <= 5 || (chr_name.size() > 5 && DEBUG_MODE)) {
			cout << "Processing chr: " << chr_name << " with " << candidate_region.size() << " regions at " << ctime(&raw_time) << endl;
			cout << endl;
		}
		for( vector< std::pair<int, int> >::iterator it = candidate_region.begin(); it != candidate_region.end(); it++ ) {
			processProperReadsBySection( samIn, samHeader, chr_name, it->first, it->second, ProperQC );
		}
		string out_name = string(output_prefix) + "proper-" + chr_name;
		printProperMapBySection( out_name.c_str() ); // out is section_name.type
		properMap.clear();
	}
	samIn.Close();
	ProperQC->PrintQCsummary();
	delete ProperQC;
	
	// process disc bam in all: discMap already initialized in generating candidate regions
	// and disc bam is also nsorted
	time_t raw_time;
	time(&raw_time);
	cout << "Processing disc reads in raw bam at: " << ctime(&raw_time) << endl;
	cout << endl;
	processDiscReads( disc_bam );
	printDiscMap( output_prefix );
}


/* inner functions ***/

/******* Proper Map related ******/

// generate chr1:100-200000 for reading bam by section
void ReadMap::setCandidateRegion( vector< std::pair<int, int> > & candidate_region, int chr_length, const char * level_name)
{
	if (LEVEL < 0) { // read all chr
		candidate_region.resize(1);
		candidate_region[0].first = 0;
		candidate_region[0].second = chr_length;
		return;
	}
	
// get info from level file
// level format: win_start	level
	ifstream level_file;
	level_file.open(level_name);
	CheckInputFileStatus(level_file, level_name);
	string line;
	int Non_level = LEVEL + 1; // if non, then what is the level?
	int prev_end = 0; // previous ending position
	int prev_st = 0; // previous start position
	while(std::getline( level_file, line )) {
		std::stringstream ss;
		ss << line;
		string st_str;
		std::getline(ss, st_str, '\t');
		string evidence_level_str;
		std::getline(ss, evidence_level_str, '\t');
		string non_level_str;
		std::getline(ss, non_level_str, '\t');
		if ( stoi(evidence_level_str) >= LEVEL || stoi(non_level_str) >= Non_level ) {
			int st = stoi(st_str);
			if ( st <= prev_end ) { // continuous, only change prev_end
				prev_end = st + win_len;
			}
			else { // not continuous, add & define new
				if ( prev_end > 0 ) {
					std::pair<int, int> new_pair;
					new_pair.first = prev_st;
					new_pair.second = prev_end;
					candidate_region.push_back(new_pair);
				}
				prev_st = st;
				prev_end = st + win_len;
			}
		}
	}
// add the last one: in case all <level followed by last coord
	if ( prev_end > 0 ) {
		std::pair<int, int> new_pair;
		new_pair.first = prev_st;
		new_pair.second = prev_end;
		candidate_region.push_back(new_pair);
	}
	
	level_file.close();
}



/**** proper related ***/
void ReadMap::processProperReadsBySection( SamFile & samIn, SamFileHeader & samHeader, std::string & chr_name, int st, int ed, QC * ProperQC)
{
	bool section_status = samIn.SetReadSection(chr_name.c_str(), st, ed);
	if (!section_status) {
		std::cerr << "ERROR: Unable to set read section: " << chr_name << ": " << st << "-" << ed << std::endl;
		exit(1);
	}

	ProperDeck pDeck(win_len, ShortClip, MinClip, MinReadLen, MinQuality, REF_SEQ);
	int Counter = 0;
	SamRecord sam_rec;
	while( samIn.ReadRecord(samHeader, sam_rec) ) {
		if (!(sam_rec.getFlag() & 0x2))
			continue;
		bool pass_qc = ProperQC->PassQC(sam_rec);
		if (!pass_qc)
			continue;
		if ( sam_rec.getInsertSize() > 0 ) // left in pair
			pDeck.Add(sam_rec);
		else {
			RetrievedIndex rv = pDeck.RetrieveIndex( sam_rec );
			addRetrievedIndexToProperMap(rv);
		}
		Counter++;
		if (Counter % 1000000 == 0)
			cout << "Proper Read by section: processed " << Counter/1000000 << " million reads!" << endl;
	}
	if (DEBUG_MODE) {
		cout << "Processed " << Counter << " reads in section: " << chr_name << ": " << st << "-" << ed << endl;
		cout << endl;
	}
}


// after retrieve info 
void ReadMap::addRetrievedIndexToProperMap( RetrievedIndex & rIndex )
{
	if (!rIndex.valid)
		return; // do not add invalid pair
		
	int first_index = (rIndex.max_position - win_len) / step_len + 1;
	if (first_index < 0) // for some positions that are at the very beginning of the chr (chr17, sample-1612, eg)
		first_index = 0;
	vector<ProperMapCell>::iterator proper_it = properMap.begin();
	proper_it += first_index;
	int walk_len = rIndex.min_position / step_len - first_index;
	
	vector<int> add_index;
	_setProperMapAddIndex( add_index, rIndex );
	
	for(int i=0; i<= walk_len; i++, proper_it++) {
	 // initialize it first: only do proper-part
		if (proper_it->size() == 0)
			proper_it->resize(4);
	  // add
	  	int n_add = 0;
	  	for(ProperMapCell::iterator sub_it = proper_it->begin(); sub_it != proper_it->end(); sub_it++, n_add++) {
	  		if ( add_index[n_add] == -1 )
	  			continue; // do not add -1
	  		if ( sub_it->size() == 0 ) // initialize
	  			n_add == 0 ? sub_it->resize(2,0) : sub_it->resize(8,0);
	  		(*sub_it)[ add_index[n_add] ]++;
	  	}
	}
}


void ReadMap::printProperMapBySection( const char* out_prefix )
{
// open out file. *.0 is proper
	vector<ofstream*> fvec;
	fvec.resize(4);
	for(int i=0; i<4; i++) {
		string out_name = string(out_prefix) + "." + std::to_string(i);
		fvec[i] = new ofstream;
		fvec[i]->open(out_name.c_str());
		CheckOutFileStatus(*fvec[i], out_name.c_str());
	}

// let's make them homogeneous so no need to care about diff index in reading them
// if any vector exists, other clip vectors exist too	
	int dist = 0;
	for(std::vector<ProperMapCell>::iterator it_print = properMap.begin(); it_print != properMap.end(); it_print++, dist++) {
		if (it_print->size() == 0) // nothing
				continue;
		int fl = 0;
		for(ProperMapCell::iterator sub_it = it_print->begin(); sub_it != it_print->end(); sub_it++, fl++) {
			*fvec[fl] << dist;
			if (sub_it->size() == 0) { // this cell does not exist but other exists
				int max_size = ( sub_it == it_print->begin() ) ? 2 : 8;
				for( int i = 0; i < max_size; i++ )
					*fvec[fl] << "\t0";
			}
			else {
				for ( vector<int>::iterator vit = sub_it->begin(); vit != sub_it->end(); vit++)
					*fvec[fl] << '\t' << *vit;
			}
			*fvec[fl] << endl;
		}
	}
// close files
	for(int i=0; i<4; i++)
		fvec[i]->close();
}



/**** ***********************utility functions ****/

// base on pairType, set add_index to proper map.
// use in RM::addRIndextoPMap
// if not suitable for add, note as -1
void ReadMap::_setProperMapAddIndex (vector<int> & add_index, RetrievedIndex & rIndex)
{
	
	add_index.resize(4);
// only proper	
	if (rIndex.type <= 2) {
		add_index[0] = (rIndex.type == 0) ? 0 : 1;
		for(int i=1; i<=3; i++)
			add_index[i] = -1;
		return;
	}
// exist clip
	add_index[0] = -1;
	char init_i = 4;
	for(int i=1; i<=3; i++) {
		init_i *= 2;
		bool mei = rIndex.type & init_i;
		bool isBegin = rIndex.type & 4;
		int actual_index = ((mei << 1) | isBegin) & 3;
		if ( rIndex.mei_on_1st )
			actual_index += 4;
		add_index[i] = actual_index;
	}
}



/******* Disc Map related ******/

// set "valid" tag based on discCandidates. (in candidate, valid = 1)
// discCandidates will be cleared after this
void ReadMap::initializeDiscMap( SamFileHeader & samHeader )
{
// clear first
	discMap.clear();

	if (discCandidates.size() == 0) {
		cerr << "ERROR: empty discCandidates!" << endl;
		exit(1);
	}

	map<string, vector< std::pair<int, int> > >::iterator discCan_it = discCandidates.begin();
	for(; discCan_it != discCandidates.end(); discCan_it++) {
	 // additional sanity check
		if (discCan_it->second.size() == 0) // no candidate on this chr
			continue;
	 // set block size
		int chr_len = std::stoi( samHeader.getSQTagValue("LN", discCan_it->first.c_str()) );
		int block_size = chr_len / step_len + 1;
		discMap[discCan_it->first].resize( block_size );
		
	// add to discMap by checking if fall within
		vector< std::pair<int, int> >::iterator sub_it = discCan_it->second.begin();
		vector<DiscMapCell>::iterator map_it = discMap[discCan_it->first].begin();
		bool reach_end = 0;
		for( int dist=0; dist < block_size; dist++, map_it++ ) {
			int cell_end = dist * step_len + win_len;
			if (!reach_end)	{
				while ( cell_end > sub_it->second ) {
					sub_it++;
					if (sub_it == discCan_it->second.end() ) {
						reach_end = 1;
						sub_it--;
						break;
					}
				}
			}
			map_it->valid = cell_end <= sub_it->second ? 1 : 0;
		}
			
		discCan_it->second.clear();
	}
	discCandidates.clear();
}


// check if a single read is in candidate: chr loc (loc + avr read length)
bool ReadMap::isSingleReadInCandidate( const char * chr_name, const int loc )
{
	string chr_str = string(chr_name);
// nothing in this chr, false
	if (discMap.find(chr_str) == discMap.end())
		return 0;
	if ( discMap[chr_str].size() == 0 )
		return 0;
// check by iterating
	vector<DiscMapCell>::iterator vit = discMap[chr_str].begin();
	int st = (loc + AvrReadLen - win_len) / step_len + 1;
	vit += st;
	int walk_len = loc / step_len - st;
	for (int i=0; i<= walk_len; i++, vit++)
		if (vit->valid)
			return 1;
// not found any available cell
	return 0;
}


// add DiscPair returned loci to discMap
void ReadMap::addLociToDiscMap( Loci & loci )
{
// decide type	
	vector<int> vec_mei;
	vec_mei.resize(3);
	for (int which_mei = 0; which_mei <3; which_mei++) {
		vec_mei[ which_mei ] = getIndexInDiscMapCellFromLoci( loci, which_mei );
	}
	
// add to map
	vector<DiscMapCell>::iterator vit = discMap[loci.chr].begin();
	int st_index = (loci.ed - win_len) / step_len + 1;
	vit += st_index;
	int walk_len = loci.st / step_len - st_index;	
	for(int dist=0; dist<= walk_len; dist++, vit++) {
		if (!vit->valid) // skip invalid windows
			continue;
		if (vit->stats.size() != 3)
			vit->stats.resize(3);
		for(int i=0; i<3; i++) {		
			if (vit->stats[i].size() != 8)
				vit->stats[i].resize(8,0);
			vit->stats[i].at( vec_mei[i] )++;
		}
	}
}


// get index from loci em + mateMap + strand
int ReadMap::getIndexInDiscMapCellFromLoci( Loci & loci , int which_mei)
{
	int trueIndex = 0;
	
	trueIndex = ((loci.type >> which_mei) & 1 );
	if (loci.sense_strand)
		trueIndex += 2;
	if (!loci.mateMap)
		trueIndex += 4;
	
	return trueIndex;
}


// process disc read (all, since it's name sorted)
void ReadMap::processDiscReads( const char * disc_bam )
{
	SamFile samIn;
	samIn.OpenForRead(disc_bam);
	if ( !samIn.IsOpen() ) {
		std::cerr << "ERROR: Unable to open disc bam: " << disc_bam << std::endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool disc_header_status = samIn.ReadHeader(samHeader);
	if (!disc_header_status) {
		cerr << "ERROR: Fail to read header: " << disc_bam << endl;
		exit(1);
	}
	
	initializeDiscMap( samHeader );

// let's process....
	DiscPair * dp = NULL;
	int Counter = 1;
	SamRecord sam_rec;
	while( samIn.ReadRecord(samHeader, sam_rec) ) {
	  // first check if exist the same
		if (dp) {
			bool same_pair = dp->IsSamePair(sam_rec);
			if (!same_pair) {
				delete dp;
				dp = NULL;
			}
		}
		if (!dp) { // 1st read, create new
			bool can1;
			bool can2;
			if ( sam_rec.getFlag() & 0x4 )
				can1 = 0;
			else
				can1 = isSingleReadInCandidate( sam_rec.getReferenceName(), sam_rec.get1BasedPosition() );
			if ( sam_rec.getFlag() & 0x8 )
				can2 = 0;
			else
				can2 = isSingleReadInCandidate( sam_rec.getMateReferenceName(), sam_rec.get1BasedMatePosition() );
			dp = new DiscPair( can1, can2, sam_rec, REF_SEQ );
		}
		else { // 2nd read
			dp->AddSecondToPair(sam_rec);
		  // add
			if ( dp->FirstValid() ) {
				Loci first_loci = dp->GetFirstLoci();
				addLociToDiscMap( first_loci );
			}
			if ( dp->SecondValid() ) {
				Loci second_loci = dp->GetSecondLoci();
				addLociToDiscMap( second_loci );
			}
		}
		Counter++;
	}
	samIn.Close();
	cout << "Processed " << Counter << " reads in disc bam: " << disc_bam << endl;
}

// print disc result out: cord-start stats...
void ReadMap::printDiscMap( const char * out_prefix )
{
// print
	for( map<string, vector<DiscMapCell> >::iterator discmap_it = discMap.begin(); discmap_it != discMap.end(); discmap_it++ ) {
		if (discmap_it->second.size() == 0) // skip no info part
			continue;	
	// initialize out file	
		vector<ofstream*> discOut;
		discOut.resize(3);
		for(int i=0; i<3; i++) {
			string fname = string(out_prefix) + "disc-" + discmap_it->first + "." + std::to_string(i);
			discOut[i] = new ofstream;
			discOut[i]->open(fname.c_str());
			CheckOutFileStatus( *discOut[int(i)], fname.c_str() );
		}
		vector<DiscMapCell>::iterator disc_it = discmap_it->second.begin();
		for( unsigned int dist = 0; dist < discmap_it->second.size(); dist++, disc_it++ ) {
			if (!disc_it->valid) // skip non-rec
				continue;
			if (disc_it->stats.size() != 3) // rarely happens when no available reads in valid region
				continue;
			for(int i=0; i<3; i++) {
				if (disc_it->stats[i].size() != 8) // skip non-existing
					continue;
				*discOut[i] << dist;
				for(vector<int>::iterator vit = disc_it->stats[i].begin(); vit != disc_it->stats[i].end(); vit++)
					*discOut[i] << "\t" << *vit;
				*discOut[i] << endl;
			}
		}
		
	// close them
		for(int i=0; i<3; i++)
			discOut[i]->close();
	}
}













































