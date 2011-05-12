/**************************************************************
 *** RHmm version 1.5.0
 ***                                                         
 *** File: cDiscrete.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
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
                cCyclicVector<cDMatrix>        mProbaMatVector;
        public :
                cDiscrete(uint theNClass, uint theNProba) ;
                virtual ~cDiscrete() ;
                virtual void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)  ;
                virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL)  ;
                virtual void InitParameters(cBaumWelchInParam &theInParam) ;
                virtual void Print() ;
                uint GetNProba() ;
                virtual void GetParam(uint theDeb, cDVector& theParam) ;
                virtual void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetNParam(void){ return mProbaMatVector[0].mNCol - 1 ; }
                void CopyDistr(cDistribution *theSrc) ;
} ;

#endif //_CDISCRETE_H_
