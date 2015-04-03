#ifndef BAMUTILITY_H
#define BAMUTILITY_H

#include <string>
#include "SamFileHeader.h"

using std::string;

// check if chr exist in bam by checking header
bool ExistChrInBam( SamFileHeader & samHeader, const char * focus_chr );

// sort bam by name, return name: (base)disc_name-nsort.bam
string SortBamByName( const char* disc_name );

int GetAvrReadLenFromBam( const char* bam );

int GetAvrInsSizeFromBam( const char* bam );

#endif

