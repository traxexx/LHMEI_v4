#include "QC.h"
#include <iostream>
#include "Globals.h"

// constructor
QC::QC( int avr_read_len, int avr_ins_size ):
	MinReadLen( avr_read_len / 2 ),
  	MinDiscIns( avr_ins_size * 3 ),
	Supplementary_Info(0),
	PCR_Duplicates(0),
	QC_Fail(0),
	Secondary_Alignment(0)
{}

// destructor
QC::~QC( void ) {}


// check if a sam record pass qc
bool QC::PassQC( SamRecord & sam_rec )
{
	int flag = sam_rec.getFlag();
	if (flag & 0x800) {
		Supplementary_Info++; // skip supplementary info
		return 0;
	}
	if (flag & 0x400) {
		PCR_Duplicates++; // skip PCR duplicates;
		return 0;
	}
	if (flag & 0x200) {
		QC_Fail++; // skip reads with QC fail;
		return 0;
	}
	if (flag & 0x100) {
		Secondary_Alignment++; // skip secondary alignment
		return 0;
	}
	return 1;
}

using std::cout;
using std::endl;
using std::cerr;
// print to std cout
void QC::PrintQCsummary()
{
	cout << "QC Summary: " << endl;
	cout << "  Supplementary_Info = " << Supplementary_Info << endl;
	cout << "  PCR_Duplicates = " << PCR_Duplicates << endl;
	cout << "  QC_Fail = " << QC_Fail << endl;
	cout << "  Secondary_Alignment = " << Secondary_Alignment << endl;
	cout << endl;
}


// check if a single record could be printed to discSam (print = 1)
bool QC::DiscSamPass( SamRecord & sam_rec )
{
// check length
	if ( sam_rec.getReadLength() < MinReadLen )
		return 0;
// check insert size
	const char * equal = sam_rec.getMateReferenceNameOrEqual();
	if ( !equal ) {
		cerr << "Warning: Unable to get mate reference name at: " << sam_rec.getReadName() << endl;
		return 0;
	}
	if ( *equal == '=' ) {
		if ( sam_rec.getInsertSize() < MinDiscIns )
			return 0;
	}
// remove QO
	simplifySamRec(sam_rec);			
// indicate print to discSam	
	return 1;
}


// OQ
void simplifySamRec( SamRecord & sam_rec )
{
	bool status = sam_rec.rmTag("OQ", 'Z');
	if (!status && DEBUG_MODE)
		std::cerr << "Warning: OQ tag not removed in " << sam_rec.getReadName() << std::endl;
	status = sam_rec.rmTag("XA", 'Z');
	if (!status && DEBUG_MODE)
		std::cerr << "Warning: XA tag not removed in " << sam_rec.getReadName() << std::endl;
	status = sam_rec.rmTag("MD", 'Z');
	if (!status && DEBUG_MODE)
		std::cerr << "Warning: MD tag not removed in " << sam_rec.getReadName() << std::endl;
	status = sam_rec.rmTag("RG", 'Z');
	if (!status && DEBUG_MODE)
		std::cerr << "Warning: RG tag not removed in " << sam_rec.getReadName() << std::endl;			
}
