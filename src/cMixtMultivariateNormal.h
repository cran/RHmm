/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cMixtMultivariateNormal.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CMIXTMULTIVARIATENORMAL_H_
#define _CMIXTMULTIVARIATENORMAL_H_
#pragma once
#include "cDistribution.h"

class cMixtMultivariateNormal : public cDistribution
{       private :
                uint    mvNClass        ;
                uint    mvNMixt         ;
                uint    mvDimObs        ;
        public :
                cDVector**     mMean   ;
                cDMatrix**     mCov    ;
                cDVector*      mp              ;
        public :
                cMixtMultivariateNormal(uint theNClass = 0, uint theNMixt = 1, uint theDimObs=1) ;
                virtual ~cMixtMultivariateNormal() ;
                void ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)  ;
                void UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba=NULL)  ;
                void InitParameters(cBaumWelchInParam &theInParam) ;
                void Print() ;
                void CopyDistr(cDistribution *theSrc) ;
                void GetParam(uint theDeb, cDVector& theParam) ;
                void SetParam(uint theDeb, cDVector& theParam) ;
                uint GetNParam(void){ return mvNMixt* mvDimObs + mvNMixt*mvDimObs*(mvDimObs+1)/2 + mvNMixt - 1 ; } ;
                uint GetNFreeParam(void){ return mvNMixt* mvDimObs + mvNMixt*mvDimObs*(mvDimObs+1)/2 + mvNMixt - 1 ; } ;
} ;
#endif //_CMIXTMULTIVARIATENORMAL_H_
