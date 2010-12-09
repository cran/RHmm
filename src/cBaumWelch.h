/**************************************************************
 *** RHmm version 1.4.4                                     
 ***                                                         
 *** File: cBaumWelch.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/09                                     
 ***                                                         
 **************************************************************/

#ifndef _CBAUMWELCH_H_
#define _CBAUMWELCH_H_
#pragma once

#include "cInParam.h"
#include "cHmm.h"

class cBaumWelch
{       private :
                        uint    mvNSample               ;
                        uint*   mvT                             ;
        public :                                                
                        cOTMatrix*      mAlpha          ;
                        cOTMatrix*      mBeta           ;
                        cOTVector*      mRho            ;
                        cOTMatrix*      mGamma          ;
                        cOTMatrix**     mXsi            ;
                        cOTMatrix*      mSumXsi         ;
                        cOTVector       mLogVrais       ;
        public :
                cBaumWelch(uint theNSample, uint* theT, uint theNClass) ;
                cBaumWelch(const cInParam &theInParam) ;
                void ForwardBackward(cOTMatrix* theCondProba, cHmm& theHMM) ;
                uint GetSampleSize(uint theN){ return mvT[theN] ;}
                virtual ~cBaumWelch() ;
} ;


#endif //_CBAUMWELCH_H_
