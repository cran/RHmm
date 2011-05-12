/**************************************************************
 *** RHmm version 1.5.0
 ***                                                         
 *** File: cMultivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CMULTIVARIATENORMAL_H_
#define _CMULTIVARIATENORMAL_H_
#pragma once

#include "cDistribution.h"
#include "SamplesUtil.h"

class cMultivariateNormal : public cDistribution
{       private :
                uint    mvNClass        ;
        public :
                cDVector*      mMean   ;
                cDMatrix*      mCov    ;
        public :
                cMultivariateNormal(uint theNClass = 0, uint theDimObs = 1) ;
                virtual ~cMultivariateNormal() ;
                virtual void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba) ;
                virtual void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL) ;
                virtual void InitParameters(cBaumWelchInParam &theInParam) ;
                virtual void Print() ;
                virtual void GetParam(uint theDeb, cDVector& theParam) ;
                virtual void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetDimObs() ;
                void CopyDistr(cDistribution* theSrc) ;
                uint GetNParam(void){ return mMean[0].mSize + mMean[0].mSize * mMean[0].mSize ; }

} ;
#endif //_CMULTIVARIATENORMAL_H_
