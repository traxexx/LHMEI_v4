#ifndef GLOBALS_H
#define GLOBALS_H

extern bool DEBUG_MODE; // if set 1 print detailed message in functions
extern bool SINGLE_SIDE; // if include single-end results in printing vcf
extern bool PSEUDO_CHR; // if 1, include pseudo chromosomes

extern int WIN;  // win length
extern int STEP; // step length
extern int LEVEL; // minimum level
extern int NON_OFFSET; // minimum # non-evidence reads = LEVEL + NON_OFFSET

#endif