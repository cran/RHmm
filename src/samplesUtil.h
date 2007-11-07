/***********************************************************
 * RHmm version 0.9.3                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/11/07                                        *
 *                                                         *
 ***********************************************************/
#ifndef _SAMPLESUTIL_H_
#define _SAMPLESUTIL_H_
#include "OTMathUtil.h"

void flatSamples(cOTVector* theInVect, uint theNSample, uint theDimObs, uint theNObsAllSamples, cOTVector& theOutVect) ;
void listSamples(cOTVector& theInVect, uint theNSample, uint theDimObs, uint* theNObsSample, cOTVector* theOutVect) ;
#endif // _SAMPLESUTIL_H_
