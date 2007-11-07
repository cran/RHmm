/***********************************************************
 * RHmm version 0.9.3                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/11/07                                        *
 *                                                         *
 ***********************************************************/
/***********************************************************
 * RHmm version 0.9.2                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/10/25                                        *
 *                                                         *
 ***********************************************************/
#ifndef _CBAUMWELCH_H_
#define _CBAUMWELCH_H_
#include "cInParam.h"
#include "cHmm.h"

class cBaumWelch
{	private :
			uint	mvNSample		;
			uint*	mvT				;
	public :						;
			cOTMatrix*	mAlpha		;
			cOTMatrix*	mBeta		;
			cOTVector*	mRho		;
			cOTMatrix*	mGamma		;
			cOTMatrix*	mXsi		;
			cOTVector	mLogVrais	;
	public :
		cBaumWelch(uint theNSample, uint* theT, uint theNClass) ;
		cBaumWelch(const cInParam &theInParam) ;
		void ForwardBackward(cOTMatrix* theCondProba, cHmm& theHMM) ;
		uint GetSampleSize(uint theN){ return mvT[theN] ;}
		virtual ~cBaumWelch() ;
} ;

#endif // _CBAUMWELCH_H_
