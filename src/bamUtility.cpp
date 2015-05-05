#include "bamUtility.h"
#include "StringBasics.h"
#include "Utilities.h"
#include "SamFile.h"
#include "Globals.h"
#include <algorithm>    // std::min_element, std::max_element

bool ExistChrInBam( SamFileHeader & samHeader, const char * focus_chr )
{
	bool focus_exist = 0;
	string focus_chr_str = string(focus_chr);
	for(int i=0; i< samHeader.getNumSQs(); i++) {  
		String chr_name = samHeader.getReferenceLabel(i);
		if ( focus_chr_str.compare(chr_name.c_str()) == 0 ) {
			focus_exist = 1;
			break;
		}
	}
	return focus_exist;
}


// sort bam by name, return name: (base)disc_name-nsort.bam
string SortBamByName( const char* disc_name )
{
	string disc_str = string(disc_name);
	string nsort_name = disc_str.substr( 0, disc_str.find_first_of('.') ) + "-nsort";
	string cmd = string("samtools sort -n ") + disc_name + " " + nsort_name;
	ExecuteCmd(cmd);
	nsort_name += ".bam";
	
	return nsort_name;
}

// set avr read len based on 1st 201~300st reads
int GetAvrReadLenFromBam( const char* bam )
{
		
// get info from sam
	SamFile samIn;
	samIn.OpenForRead(bam);
	if ( !samIn.IsOpen() ) {
		std::cerr << "ERROR: Unable to open AvrReadLen bam: " << bam << std::endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool bam_header_status = samIn.ReadHeader(samHeader);
	if (!bam_header_status) {
		std::cerr << "ERROR: Fail to read header: " << bam << std::endl;
		exit(1);
	}
	
// let's process....
	int Counter = 0;
	SamRecord sam_rec;
	int avr_read_len = 0;
	while( samIn.ReadRecord(samHeader, sam_rec) ) {
		Counter++;
		if (Counter > 300)
			break;
		else if (Counter < 200)
			continue;
		avr_read_len += sam_rec.getReadLength();
	}
	samIn.Close();

// get avr, return	
	avr_read_len /= (Counter - 200);	
	return avr_read_len;
}


// avr read pair insert size: based on 1st 201~300 0x2 reads
int GetAvrInsSizeFromBam( const char* bam )
{
	SamFile samIn;
	samIn.OpenForRead(bam);
	if ( !samIn.IsOpen() ) {
		std::cerr << "ERROR: Unable to open AvrReadLen bam: " << bam << std::endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool bam_header_status = samIn.ReadHeader(samHeader);
	if (!bam_header_status) {
		std::cerr << "ERROR: Fail to read header: " << bam << std::endl;
		exit(1);
	}
	
// process
	int Counter = 0;
	SamRecord sam_rec;
	int avr_ins_size = 0;
	while( samIn.ReadRecord(samHeader, sam_rec) ) {
		if ( !(sam_rec.getFlag() & 0x2) ) // skip disc
			continue;
		Counter++;
		if (Counter > 300)
			break;
		else if (Counter < 200)
			continue;
		avr_ins_size += abs( sam_rec.getInsertSize() );
	}
	samIn.Close();

// get avr, return	
	avr_ins_size /= (Counter - 200);	
	return avr_ins_size;
}

// get 100 regions in chr20 then calculate depth
float EstimateBamDepth( string & bam, int & readlen)
{
// parameters
	int pad = 200000;
	int cwin = 2000;
	int nwin = 60;
	int minlen = cwin * nwin;

// bam sanity check
	SamFile samIn;
	samIn.OpenForRead( bam.c_str() );
	SamFileHeader samHeader;
	bool bam_header_status = samIn.ReadHeader(samHeader);
	if (!bam_header_status) {
		std::cerr << "ERROR: Fail to read header: " << bam << std::endl;
		exit(1);
	}
	string bai_name = bam + ".bai";
	bool bai_status = samIn.ReadBamIndex( bai_name.c_str() );
	if (!bai_status) {
		std::cerr << "ERROR: Invalid bam index: " << bai_name << std::endl;
		exit(1);
	}

// get ref chr length & set regions
	int reflen = atoi( samHeader.getSQTagValue("LN", REF_CHR.c_str()) );
	int trimlen = reflen - pad * 2;
	if ( trimlen < minlen ) {
		std::cerr << "ERROR: ref-chr: " << REF_CHR << " is too short to estimate read depth. Are you using the right ref chromosome?" << std::endl;
		exit(1);
	}
	
// count read number
	int times = trimlen / minlen;
	vector<int> vlc;
	vlc.resize( nwin );
	for( int seg = 0; seg < nwin; seg++ ) {
		bool section_status = samIn.SetReadSection(REF_CHR.c_str(), pad + seg*times, pad + seg*times + cwin);
		if (!section_status) {
			std::cerr << "ERROR: Unable to set read section: " << REF_CHR << ": " << pad + seg*times << "-" << pad + seg*times + cwin << std::endl;
			exit(1);
		}
		int lc = 0;
		SamRecord sam_rec;
		while( samIn.ReadRecord(samHeader, sam_rec) ) {
			lc++;
		}
		vlc[ seg ] = lc;
	}
	samIn.Close();
	int rcount = 0;
	for( vector<int>::iterator it = vlc.begin(); it != vlc.end(); it++ ) {
		rcount += *it;
	}
	// remove max & min
	rcount -= ( *std::min_element( vlc.begin(), vlc.end() ) + *std::max_element( vlc.begin(), vlc.end() ) );
	
// depth
	float depth;
	int region_size = minlen - 2 * cwin - 2 * readlen * nwin;
	depth = float(rcount) / (region_size / WIN) / 2;
	return depth;
}

