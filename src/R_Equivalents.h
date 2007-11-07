/***********************************************************
 * RHmm version 0.9.3                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/11/07                                        *
 *                                                         *
 ***********************************************************/
#ifndef _REQUIVALENTS_
#define _REQUIVALENTS_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Rprintf	printf
inline void GetRNGstate(void) {srand( (unsigned)time( NULL ) ); }
inline double unif_rand(void) { return (double)rand()/(double)RAND_MAX ; }
inline void PutRNGstate(void){}
#endif //_REQUIVALENTS_
