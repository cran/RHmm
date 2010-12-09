/**************************************************************
 *** RHmm version 1.4.4                                     
 ***                                                         
 *** File: cViterbi.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/09                                     
 ***                                                         
 **************************************************************/

#ifndef _CVITERBI_H_
#define _CVITERBI_H_
#pragma once
#include "cInParam.h"
#include "cHmm.h"
#include "cDistribution.h"

class cViterbi
{       public :
                uint            **mSeq          ;
                cOTVector       mLogProb        ;
        public :
                cViterbi(cInParam &theInParam) ;
                ~cViterbi() ;
                void ViterbiPath(cInParam& theInParam, cHmm& theHMM) ;
} ;
#endif //_CVITERBI_H_
