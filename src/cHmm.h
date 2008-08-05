/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
 *                                                         *
 ***********************************************************/
#ifndef _CHMM_H_
#define _CHMM_H_
#include "hmm.h"
#include "cinparam.h"

class cDistribution ;
class cHmm
{	public :
		distrDefinitionEnum	mDistrType	;
		cOTVector			mInitProba	;
		cOTMatrix			mTransMat	;	
		cDistribution*		mDistrParam	;
	public :
		cHmm(distrDefinitionEnum theDistrType, uint theNClass , uint theDimObs=1 , uint theNMixt=0, uint theNProba=0) ;
		cHmm(const cInParam &theInParam) ;
		virtual ~cHmm() ;
		cHmm & operator =(cHmm& theSrc) ;
		void CopyHmm(cHmm& theSrc) { *this = theSrc ; } ;
		void Print() ;
		uint GetNParam(void) ;
		void SetParam(cOTVector& theParam) ;
		void GetParam(cOTVector& theParam) ;		
} ;



#endif _CHMM_H_//