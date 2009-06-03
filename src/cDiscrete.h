/**************************************************************
 *** RHmm version 1.3.0                                      
 ***                                                         
 *** File: cDiscrete.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2009/06/03                                      
 ***                                                         
 **************************************************************/

#ifndef _CDISCRETE_H_
#define _CDISCRETE_H_
#include "cDistribution.h"
class cDiscrete : public cDistribution 
{	private :
		uint		mvNClass	;
	public :
		cOTVector*	mProba		;
	public :
		cDiscrete(uint theNClass, uint theNProba) ;
		virtual ~cDiscrete() ;
		virtual void ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba)  ;
		virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba=NULL)  ;
		virtual void InitParameters(cBaumWelchInParam &theInParam) ;
		virtual void Print() ;
		uint GetNProba() ;
		virtual void GetParam(uint theDeb, cOTVector& theParam) ;
		virtual void SetParam(uint theDeb, cOTVector& theParam) ;
		uint GetNParam(void){ return mProba[0].mSize - 1 ; }
		void CopyDistr(cDistribution *theSrc) ;
} ;

#endif //_CDISCRETE_H_
