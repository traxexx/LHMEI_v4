#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "SamFile.h"
#include "SamFileHeader.h"
#include "MeiCoord.h"
#include "QC.h"

// only take single-chr out: use in single-chr LHMEI
void TakeChrOut( const char* inBam, const char* outSam, const char* ctrlChr );

// take ref-chr out + print disc sam
void PreProcessBam( const char* inBam, const char* outSam, const char* ctrlChr,
	 const char * disc_name, const char * mei_coord_list, int avr_read_len, int avr_ins_size);

// generate level list from disc.sam
void GenerateLevelListFromDiscBam( const char * disc_name, const char * list_prefix, int step, int win );

	 
/************ inner functions **************/

// non-ref sections
void processNonRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader, QC & BamQC, MeiCoord & mei_coord, string & REF_CHR );


// ref section (print all to samOut)
void processRefSection( SamFile & samIn, SamFile & samOut, SamFile & discSam, SamFileHeader & samHeader, QC & BamQC, MeiCoord & mei_coord );

#endif
