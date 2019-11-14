/**************************************************************
 *** RHmm version 1.3.4                                      
 ***                                                         
 *** File: cViterbi.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2010/11/14                                      
 ***                                                         
 **************************************************************/

#ifndef _CVITERBI_H_
#define _CVITERBI_H_
#include "cInParam.h"
#include "cHmm.h"
#include "cDistribution.h"

class cViterbi
{	public :
		uint		**mSeq		;
		cOTVector	mLogProb	;
	public :
		cViterbi(cInParam &theInParam) ;
		~cViterbi() ;
		void ViterbiPath(cInParam& theInParam, cHmm& theHMM) ;
} ;

#endif // _CVITERBI_H_
