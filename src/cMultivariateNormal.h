/**************************************************************
 *** RHmm version 1.3.0                                      
 ***                                                         
 *** File: cMultivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2009/06/03                                      
 ***                                                         
 **************************************************************/

#ifndef _CMULTIVARIATENORMAL_H_
#define _CMULTIVARIATENORMAL_H_
#include "cDistribution.h"
#include "SamplesUtil.h"

class cMultivariateNormal : public cDistribution
{	private :
		uint	mvNClass	;
	public :
		cOTVector*	mMean	;
		cOTMatrix*	mCov	;
	public :
		cMultivariateNormal(uint theNClass = 0, uint theDimObs = 1) ;
		virtual ~cMultivariateNormal() ;
		virtual void ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba) ;
		virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba=NULL) ;
		virtual void InitParameters(cBaumWelchInParam &theInParam) ;
		virtual void Print() ;
		virtual void GetParam(uint theDeb, cOTVector& theParam) ;
		virtual void SetParam(uint theDeb, cOTVector& theParam) ;
		uint GetDimObs() ;
		void CopyDistr(cDistribution* theSrc) ;
		uint GetNParam(void){ return mMean[0].mSize + mMean[0].mSize * mMean[0].mSize ; }

} ;

#endif //_CMULTIVARIATENORMAL_H_