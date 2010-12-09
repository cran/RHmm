/**************************************************************
 *** RHmm version 1.4.4                                     
 ***                                                         
 *** File: REquivalents.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/09                                     
 ***                                                         
 **************************************************************/

#ifndef _REQUIVALENTS_H_
#define _REQUIVALENTS_H_
#pragma once
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#ifdef __SUNPRO_CC
    #define log std::log
    #define exp std::exp
    #define fabs std::fabs
    #define printf std::printf
    #define pow std::pow
    #define sqrt std::sqrt
#endif __SUNPRO_CC

#define Rprintf printf
inline void GetRNGstate(void) { std::srand( (unsigned)time( NULL ) ); }
inline double unif_rand(void) { return (double)std::rand()/(double)RAND_MAX ; }
inline void PutRNGstate(void){}

#endif //_REQUIVALENTS_H_
