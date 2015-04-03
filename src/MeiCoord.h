#ifndef MEICOORD_H
#define MEICOORD_H

#include "SamRecord.h"
#include <vector>
#include <map>
#include <string>
#include <iterator>

using std::vector;
using std::map;
using std::string;

typedef struct {
	int start;
	int end;
} Coord;

typedef vector<Coord>::iterator CordPtr;

typedef vector< vector<Coord> > MultiCoord;

typedef map< string, MultiCoord > MultiCoordMap;

class MeiCoord
{
  public:
	MeiCoord( const char * mei_coord_list );
	~MeiCoord( void );
	
	void ResetVecPtr( string & current_chr ); // do it when change chr
  	void SetEMtag( SamRecord & sam_rec );
  	
  private:
	MultiCoordMap RefCoord;
	
	vector<CordPtr> CordPtrVec; // for set EM tag of disc
	vector<CordPtr> EndPtrVec; // for end purpose
	
	bool No_rec_on_this_chr;
};

#endif