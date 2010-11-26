/**************************************************************
 *** RHmm version 1.4.2                                     
 ***                                                         
 *** File: REquivalents.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/11/26                                     
 ***                                                         
 **************************************************************/

#ifndef _REQUIVALENTS_H_
#define _REQUIVALENTS_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Rprintf printf
inline void GetRNGstate(void) {srand( (unsigned)time( NULL ) ); }
inline double unif_rand(void) { return (double)rand()/(double)RAND_MAX ; }
inline void PutRNGstate(void){}
#endif //_REQUIVALENTS_H_
