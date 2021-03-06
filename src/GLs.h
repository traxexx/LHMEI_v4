#ifndef GLS_H
#define GLS_H

#include <vector>
#include <string>

using std::vector;
using std::string;

// utility functions in operating GL

float GetGLfromCounts( vector<int> & counts, vector<float> & ref );

float SumGL( float original, float single_gl );

float MinusGL( float original, float compensate);

int GetVecSum( vector<int> & counts );

int GetNumberOfZerosInVec( vector<int> & counts );

/*** PL related *****/

float GetProbFromGLs( float BaseGL, float AddGL );

int GetVariantQuality( vector<float> & GL );

int GetAlleleDosage( vector<float> & GL );

string GetGenotype( vector<float> & GL );

void SetPLsFromGL( vector<int> & PL, vector<float> & GL );

int GetGenotypeQuality( vector<float> & GL );


#endif