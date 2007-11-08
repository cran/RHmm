/***********************************************************
 * RHmm version 0.9.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/11/08                                        *
 *                                                         *
 ***********************************************************/
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
