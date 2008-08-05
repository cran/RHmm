/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
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
