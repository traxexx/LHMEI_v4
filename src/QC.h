#ifndef QC_H
#define QC_H

#include "SamFile.h"

class QC
{
  public:	
	QC( int avr_read_len, int avr_ins_size );
	~QC( void );
	bool PassQC( SamRecord & sam_rec );
	bool DiscSamPass( SamRecord & sam_rec );
	void PrintQCsummary();

  private:
  	const int MinReadLen;
  	const int MinDiscIns;
	int Supplementary_Info;
	int PCR_Duplicates;
	int QC_Fail;
	int Secondary_Alignment;	
};


// remove QO tags for disc sam
void simplifySamRec( SamRecord & sam_rec );

#endif
