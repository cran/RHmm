/**************************************************************
 *** RHmm version 1.3.0                                      
 ***                                                         
 *** File: Kmeans.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2009/04/27                                       
 ***                                                         
 **************************************************************/

#ifndef _KMEANS_H_
#define _KMEANS_H_
#ifndef _RDLL_

#include "OTMathUtil.h"
#include "REquivalents.h"
#ifndef uint
	typedef unsigned int uint ;
#endif //uint


extern void mkmeans(cOTVector& theYt, uint theNClass, int* theSeq) ;
extern void mkmeans(cOTVector& theYt, uint theNClass, uint theDimObs, int* theSeq) ;
#endif //_RDLL_
#endif // _KMEANS_H_


