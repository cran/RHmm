/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
 *                                                         *
 ***********************************************************/
#ifndef _CUNIVARIATENORMAL_H_
#define _CUNIVARIATENORMAL_H_
#include "cdistribution.h"

class cUnivariateNormal : public cDistribution
{	public :
		cOTVector	mMean	;
		cOTVector	mVar	;
	public :
		cUnivariateNormal(uint theNClass = 0) ;
		virtual ~cUnivariateNormal() ;
		virtual void ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba) ;
		virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba=NULL) ;
		virtual void InitParameters(cBaumWelchInParam &theInParam) ;
		virtual void Print() ;
		virtual void GetParam(uint theDeb, cOTVector& theParam) ;
		virtual void SetParam(uint theDeb, cOTVector& theParam) ;
		virtual uint GetNParam(void){ return 2 ; }
		void CopyDistr(cDistribution* theSrc) ;
} ;

#endif //_CUNIVARIATENORMAL_H_
