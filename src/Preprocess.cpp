#include "StringArray.h"
#include "StringHash.h"
#include "Parameters.h"
#include "Error.h"

#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <fstream>

#include "Preprocess.h"
#include "QC.h"
#include "Utilities.h"
#include "Globals.h"

using std::string;
using std::ifstream;
using std::ofstream;

// check if out/in bam is open. also check .bai status
void SanityCheckBams( SamFile & samIn, SamFile & samOut, bool & bai_status )
{
	if ( !samIn.IsOpen() ) {
		std::cerr << "ERROR: Unable to open input bam " << std::endl;
		exit(1);
	}
	
	if ( !samOut.IsOpen() ) {
		std::cerr << "ERROR: Unable to open output bam " << std::endl;
		exit(1);
	}
	
	if ( !bai_status ) {
		std::cerr << "ERROR: Unable to read bam index, is it indexed?" << std::endl;
		exit(1);
	}
}


void TakeChrOut( const char* inBam, const char* outSam, const char* ctrlChr )
{
	SamFile samIn, samOut;
	samIn.OpenForRead(inBam);
	samOut.OpenForWrite(outSam);
	
	std::string bai = std::string(inBam) + ".bai";
	bool bai_status = samIn.ReadBamIndex( bai.c_str() );
	
	SanityCheckBams( samIn, samOut, bai_status );
	
	SamFileHeader samHeader;
	bool in_bam_header_status = samIn.ReadHeader(samHeader);
	if (!in_bam_header_status) {
		std::cerr << "ERROR: Fail to read header: " << inBam << std::endl;
		exit(1);
	}
	samOut.WriteHeader(samHeader);

// list of non_ctrl chr
	std::string REF_CHR = std::string(ctrlChr);
	std::vector<std::string> chrNameVec;
	for(int i=0; i< samHeader.getNumSQs(); i++) {  
		String val = samHeader.getReferenceLabel(i);
		std::string str_val = std::string(val.c_str());
		chrNameVec.push_back(str_val);
	}
	
// process bam
	for(int i=0; i< samHeader.getNumSQs(); i++) {
		std::string current_chr = chrNameVec[i];	
		bool section_status = samIn.SetReadSection(current_chr.c_str());
		if (!section_status) {
			std::cerr << "ERROR: Unable to set read section to chr " << current_chr << std::endl;
			exit(1);
		}
		
		int Counter = 0;
		SamRecord samRecord;
		
	// ref section
		if ( current_chr.compare(REF_CHR) == 0) {
			 while(samIn.ReadRecord(samHeader, samRecord)) {
			 	Counter++;
         		samOut.WriteRecord(samHeader, samRecord);
			 }
		}
		else { // non-ref section
 			while(samIn.ReadRecord(samHeader, samRecord)) {
         		Counter++;
         		int flag = samRecord.getFlag();
         		if (flag & 0x2)  continue;
         		if (flag & 0x8)  continue; // mate unmap
         		std::string mate_chr = samRecord.getMateReferenceName();
         		if (mate_chr.compare(REF_CHR) == 0)
	         		samOut.WriteRecord(samHeader, samRecord);
      		}		
		}
		if (current_chr.length() <= 5)
			std::cout << "Processed chr " << current_chr << ", total = " << Counter << " reads." << std::endl;
	}

}


void PreProcessBam( const char* inBam, const char* outSam, const char* ctrlChr,
	const char * disc_name, const char * mei_coord_list, int avr_read_len, int avr_ins_size)
{
	SamFile samIn, samOut;
	samIn.OpenForRead(inBam);
	samOut.OpenForWrite(outSam);
	
	std::string bai = std::string(inBam) + ".bai";
	bool bai_status = samIn.ReadBamIndex( bai.c_str() );
	
	SanityCheckBams( samIn, samOut, bai_status );
	
	SamFileHeader samHeader;
	bool in_bam_header_status = samIn.ReadHeader(samHeader);
	if (!in_bam_header_status) {
		std::cerr << "ERROR: Fail to read header: " << inBam << std::endl;
		exit(1);
	}
	samOut.WriteHeader(samHeader);
	
	SamFile discSam;
	discSam.OpenForWrite( disc_name );
	if ( !discSam.IsOpen() ) {
		std::cerr << "Unable to open " << disc_name << std::endl;
		exit(1);
	}
	discSam.WriteHeader( samHeader);

// initialize mei-coord-ref
	MeiCoord mei_coord( mei_coord_list );

// list of non_ctrl chr
	std::string REF_CHR = std::string(ctrlChr);
	std::vector<std::string> chrNameVec;
	for(int i=0; i< samHeader.getNumSQs(); i++) {  
		String val = samHeader.getReferenceLabel(i);
		std::string str_val = std::string(val.c_str());
		chrNameVec.push_back(str_val);
	}

	QC BamQC( avr_read_len, avr_ins_size );
// process bam by section
	std::cout << std::endl;
	std::cout << "Pre-process raw bam..." << std::endl;
	for(int i=0; i< samHeader.getNumSQs(); i++) {
		std::string current_chr = chrNameVec[i];	
		bool section_status = samIn.SetReadSection(current_chr.c_str());
		if (!section_status) {
			std::cerr << "ERROR: Unable to set read section to chr " << current_chr << std::endl;
			exit(1);
		}
		
		mei_coord.ResetVecPtr( current_chr );
		if ( current_chr.compare(REF_CHR) == 0) { // ref section
			std::cout << "Working on ctrl chr: " << current_chr << std::endl;
			processRefSection( samIn, samOut, discSam, samHeader, BamQC, mei_coord );
		}
		else { // non ref
			std::cout << "Working on chr: " << current_chr << std::endl;
			processNonRefSection( samIn, samOut, discSam, samHeader, BamQC, mei_coord, REF_CHR );
		}
	}

	BamQC.PrintQCsummary();
}


/* generate level list from disc.sam
   levels:
  			0 windows with unmap/disc
  			1 windows with 1 MEI disc
  			2,3,4,5 etc
   ADD BY MATE (EM is for the current read): so not that accurate but who cares
   */
void GenerateLevelListFromDiscBam( const char * disc_name, const char * list_prefix, int step, int win )
{
	if (DEBUG_MODE)
		std::cout << "Generating level list from disc bam..." << std::endl;
	int AvrReadLen = 100; // let's temporarily set this

	SamFile samIn;
	samIn.OpenForRead(disc_name);
	if ( !samIn.IsOpen() ) {
		std::cerr << "ERROR: Can't open " << disc_name << std::endl;
		exit(1);
	}
	SamFileHeader samHeader;
	bool disc_header_status = samIn.ReadHeader(samHeader);
	if (!disc_header_status) {
		std::cerr << "ERROR: Fail to read header: " << disc_name << std::endl;
		exit(1);
	}
	
// generate list for all then print out by chr

  // generate chr map
  	map<string, vector<int> > levelMap;
  	for(int i=0; i< samHeader.getNumSQs(); i++) {
  		String chr_Str = samHeader.getReferenceLabel(i);
  		string chr_str = string(chr_Str.c_str());
  		int chr_size = std::stoi( samHeader.getSQTagValue("LN", chr_str.c_str()) );
		int win_count = chr_size / step + 1;
		levelMap[chr_str].resize(win_count, -1);
  	}
	
  // add level from bam
	SamRecord sam_rec;
	int line_count = 0;
	while(samIn.ReadRecord(samHeader, sam_rec)) {
		line_count++;
		if (sam_rec.getFlag() & 0x8) // skip unmap
			continue;
		string current_chr = sam_rec.getMateReferenceName();
		int current_index_min = (sam_rec.get1BasedMatePosition() + AvrReadLen - win) / step + 1;
		if (current_index_min < 0)
			current_index_min = 0;
		int current_index_max = sam_rec.get1BasedMatePosition() / step;
		if ( current_index_max < 0 ) {
			std::cerr << "Warning: negateive mapped position at: " << sam_rec.getReadName() << ". Skipped!" << std::endl;
			continue;	
		}
		if ( current_index_max >= int(levelMap[current_chr].size()) ) {
			std::cerr << "Warning: Mate of : " << sam_rec.getReadName() << " at: " << current_chr << "-" << sam_rec.get1BasedMatePosition() << " out of chr boundary. Skipped!" << std::endl;
			continue;
		}
		int * tag = sam_rec.getIntegerTag("EM");
		bool em_type;
		if (tag)
			em_type = *tag > 0 ? 1 : 0;
		else // no EM (that's not error, including some disc & unmap)
			em_type = 0;
		vector<int>::iterator lmap = levelMap[current_chr].begin();
		lmap += current_index_min;
		if (em_type) { // add to map
			for(int i = current_index_min; i <= current_index_max; i++, lmap++) {
				(*lmap)++;
			}
		}
		else { // can only convert -1 to 0
			for(int i = current_index_min; i <= current_index_max; i++, lmap++) {
//std::cout << line_count<< ": " << current_chr << " : " << i << ", size = " << levelMap[current_chr].size() << std::endl;
				if ( (*lmap) < 0 )
					(*lmap) = 0;
			}
		}
	}


  // print level info
	for( map<string, vector<int> >::iterator mit = levelMap.begin(); mit != levelMap.end(); mit++ ) {
		ofstream list_file;
		string list_name = string(list_prefix) + "." + mit->first;
		list_file.open(list_name.c_str());
		CheckOutFileStatus(list_file, list_name.c_str());
		vector<int>::iterator it = mit->second.begin();
		for( unsigned int dist=0; dist < mit->second.size(); dist++, it++ ) {
			if (*it > -1) {
				list_file << dist*step << "\t" << *it << std::endl;
			}
		}
		list_file.close();
	}
}



/*************** inner function ***************/

// only print: related to REF-CHR to samOut & disc to discSam
void processNonRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader,
	QC & BamQC, MeiCoord & mei_coord, string & REF_CHR )
{
	int Counter = 0;
	SamRecord sam_rec;

	while(samIn.ReadRecord(samHeader, sam_rec)) {
		Counter++;

		if ( sam_rec.getFlag() & 0x2 )
			continue;
		bool qc_pass = BamQC.PassQC( sam_rec );
		if ( !qc_pass ) 
			continue;
		if ( sam_rec.getFlag() & 0x8 ) {
			if ( sam_rec.getFlag() & 0x4 )
				continue;
			discSam.WriteRecord(samHeader, sam_rec);
			continue;
		}
		std::string mate_chr = sam_rec.getMateReferenceName();
		if (mate_chr.compare(REF_CHR) == 0)
			samOut.WriteRecord(samHeader, sam_rec);
		if ( sam_rec.getFlag() & 0x4 ) {
			discSam.WriteRecord(samHeader, sam_rec);
			continue;
		}	
	
	// reset qual if irregular chr? should we do that?
		bool disc_pass = BamQC.DiscSamPass ( sam_rec );
		if ( disc_pass ) {
			mei_coord.SetEMtag( sam_rec );
			discSam.WriteRecord(samHeader, sam_rec);
		}
	}
	std::cout << "  Total = " << Counter << " reads processed!" << std::endl;
}


// print all to samOut & disc to discSam
void processRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader,
	QC & BamQC, MeiCoord & mei_coord )
{
	int Counter = 0;
	SamRecord sam_rec;
	
	while(samIn.ReadRecord(samHeader, sam_rec)) {
		Counter++;
		if (sam_rec.getFlag() & 0x4) {
			if ( sam_rec.getFlag() & 0x8 )
				continue;
			discSam.WriteRecord(samHeader, sam_rec);
			samOut.WriteRecord(samHeader, sam_rec);
			continue;
		}
		samOut.WriteRecord(samHeader, sam_rec);
		
		if ( sam_rec.getFlag() & 0x2 ) 
			continue;
		bool qc_pass = BamQC.PassQC( sam_rec );
		if ( !qc_pass ) 
			continue;

		if ( sam_rec.getFlag() & 0x8 ) {
			discSam.WriteRecord(samHeader, sam_rec);
			continue;
		}	
	
	// reset qual if irregular chr? should we do that?
		bool disc_pass = BamQC.DiscSamPass( sam_rec );
		if ( disc_pass ) {
			mei_coord.SetEMtag( sam_rec );
			discSam.WriteRecord(samHeader, sam_rec);
		}
	}
	std::cout << "  Total = " << Counter << " reads processed!" << std::endl;
}


















