/**************************************************************
 *** RHmm version 1.5.0
 ***                                                         
 *** File: cUnivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CUNIVARIATENORMAL_H_
#define _CUNIVARIATENORMAL_H_
#pragma once
#include "cDistribution.h"

class cUnivariateNormal : public cDistribution
{       public :
                cDVector       mMean   ;
                cDVector       mVar    ;
        public :
                cUnivariateNormal(uint theNClass = 0) ;
                virtual ~cUnivariateNormal() ;
                virtual void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba) ;
                virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL) ;
                virtual void InitParameters(cBaumWelchInParam &theInParam) ;
                virtual void Print() ;
                virtual void GetParam(uint theDeb, cDVector& theParam) ;
                virtual void SetParam(uint theDeb, cDVector& theParam) ;
                virtual uint GetNParam(void){ return 2 ; }
                void CopyDistr(cDistribution* theSrc) ;
} ;
#endif //_CUNIVARIATENORMAL_H_
