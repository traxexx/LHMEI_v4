#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>

extern bool DEBUG_MODE; // if set 1 print detailed message in functions
extern bool SINGLE_SIDE; // if include single-end results in printing vcf
extern bool PSEUDO_CHR; // if 1, include pseudo chromosomes
extern bool PRINT_NON_VARIANT; // if 1, print GT=0 windows. For debug use
extern bool DEPTH_FILTER; // filter variants by depth

extern int WIN;  // win length
extern int STEP; // step length
extern int LEVEL; // minimum level
extern int NON_OFFSET; // minimum # non-evidence reads = LEVEL + NON_OFFSET

extern std::string REF_CHR;

extern float DEPTH;
//extern int READLEN;

#endif