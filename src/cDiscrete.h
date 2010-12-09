/**************************************************************
 *** RHmm version 1.4.4                                     
 ***                                                         
 *** File: cDiscrete.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/09                                     
 ***                                                         
 **************************************************************/

#ifndef _CDISCRETE_H_
#define _CDISCRETE_H_
#pragma once

#include "cDistribution.h"
class cDiscrete : public cDistribution 
{       private :
                uint            mvNClass        ;
        public :
                cOTVector*      mProba          ;
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
