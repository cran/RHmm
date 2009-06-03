/**************************************************************
 *** RHmm version 1.3.0                                      
 ***                                                         
 *** File: cMixtMultivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2009/06/03                                      
 ***                                                         
 **************************************************************/

#ifndef _CMIXTMULTIVARIATENORMAL_H_
#define _CMIXTMULTIVARIATENORMAL_H_
#include "cDistribution.h"

class cMixtMultivariateNormal : public cDistribution
{	private :
		uint	mvNClass	;
		uint	mvNMixt		;
		uint	mvDimObs	;
	public :
		cOTVector**	mMean	;
		cOTMatrix**	mCov	;
		cOTVector*	mp		;
	public :
		cMixtMultivariateNormal(uint theNClass = 0, uint theNMixt = 1, uint theDimObs=1) ;
		virtual ~cMixtMultivariateNormal() ;
		virtual void ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba)  ;
		virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba=NULL)  ;
		virtual void InitParameters(cBaumWelchInParam &theInParam) ;
		virtual void Print() ;
		virtual void CopyDistr(cDistribution *theSrc) ;
		virtual void GetParam(uint theDeb, cOTVector& theParam) ;
		virtual void SetParam(uint theDeb, cOTVector& theParam) ;
		uint GetNParam(void){ return mvNMixt* mvDimObs + mvNMixt*mvDimObs*(mvDimObs+1)/2 + mvNMixt - 1 ; } ;
} ;

#endif //_CMIXTMULTIVARIATENORMAL_H_
