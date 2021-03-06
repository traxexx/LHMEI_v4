#ifndef ANCHORRANK_H
#define ANCHORRANK_H

#include <vector>

using std::vector;

// record window merge info for print the single cell
class SingleCellPrint
{
  public:
  	SingleCellPrint( int win_index, vector<float> & GL, vector<int> & counts ); // constructor
  	~SingleCellPrint();
  	void UpdateCoordWithNew( int win_index );
  	void UpdateBothEnd( vector<int> & counts );
  	void UpdateWithEqual( int win_index, vector<int> & counts );
  	
	int wcount; // # merge window
	int central; // var position to report
	int anchor_end; // use in add
	int ci; // confidence interval
	int var_end;
	int gq_peak; // highest g qual
	bool both_end; // if both end anchor present?
};

// rule: dosage --> posterior-variant --> %(disc + clip + unmap) --> less %proper -->depth --> use anchor
//bool CompareDosage( vector<int> & anchorGL, vector<int> & newGL );

int GetVariantPosterior( vector<float> & GL);

int GetSupportReadFraction( vector<int> & counts, int depth );

int GetProperReadFraction( vector<int> & counts, int depth );

int getSumSupportClips( vector<int> & counts );

int getSumSupportDiscs( vector<int> & counts);

int getSumSupportUnmaps( vector<int> & counts);

#endif