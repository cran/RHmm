/**************************************************************
 *** RHmm version 1.5.0
 ***                                                         
 *** File: cMixtUnivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CMIXTUNIVARIATENORMAL_H_
#define _CMIXTUNIVARIATENORMAL_H_
#pragma once

#include "cDistribution.h"

class cMixtUnivariateNormal : public cDistribution
{       private :
                uint    mvNClass        ;
                uint    mvNMixt         ;
        public :
                cDVector*      mMean   ;
                cDVector*      mVar    ;
                cDVector*      mp              ;
        public :
                cMixtUnivariateNormal(uint theNClass = 0, uint theNMixt = 1) ;
                virtual ~cMixtUnivariateNormal() ;
                virtual void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)  ;
                virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL)  ;
                virtual void InitParameters(cBaumWelchInParam &theInParam) ;
                virtual void Print() ;
                virtual void CopyDistr(cDistribution *theSrc) ;
                virtual void GetParam(uint theDeb, cDVector& theParam) ;
                virtual void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetNParam(void){ return mvNMixt * 3 - 1  ; }
} ;
#endif //_CMIXTUNIVARIATENORMAL_H_
