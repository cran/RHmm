/**************************************************************
 *** RHmm version 1.4.2                                     
 ***                                                         
 *** File: SamplesUtil.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/11/26                                     
 ***                                                         
 **************************************************************/

#ifndef _SamplesUtil_H_
#define _SamplesUtil_H_
#include "OTMathUtil.h"

void flatSamples(cOTVector* theInVect, uint theNSample, uint theDimObs, uint theNObsAllSamples, cOTVector& theOutVect) ;
void listSamples(cOTVector& theInVect, uint theNSample, uint theDimObs, uint* theNObsSample, cOTVector* theOutVect) ;
#endif // _SamplesUtil_H_
