/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
 *                                                         *
 ***********************************************************/
#ifndef _CDISTRIBUTION_H_
#define _CDISTRIBUTION_H_
#include "cbaumwelchinparam.h"
#include "cbaumwelch.h"


class cDistribution  
{	public :
		virtual void ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba)=0 ;
		virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba=NULL)=0 ;
		virtual void InitParameters(cBaumWelchInParam& theInParam)=0 ;
		virtual void Print() = 0 ;
		virtual void CopyDistr(cDistribution* theSrc)=0 ;
		virtual uint GetNParam(void)=0;
		virtual void GetParam(uint theDeb, cOTVector& theParam)=0 ;
		virtual void SetParam(uint theDeb, cOTVector& theParam)=0 ;
#ifndef _RDLL_
		void KMeans(cOTVector& theYt, uint theNClass, int* theSeq) {
				mkmeans(theYt, theNClass, theSeq) ;
		} ;
		void KMeans(cOTVector& theYt, uint theNClass, uint theDimObs, int* theSeq) {
				mkmeans(theYt, theNClass, theDimObs, theSeq) ;
		} ;
#endif // _RDLL_

} ;

#endif //_CDISTRIBUTION_H_