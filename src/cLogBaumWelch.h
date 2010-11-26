/**************************************************************
 *** RHmm version 1.4.2                                     
 ***                                                         
 *** File: cLogBaumWelch.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/11/26                                     
 ***                                                         
 **************************************************************/

#ifndef _CLOGBAUMWELCH_H_
#define _CLOGBAUMWELCH_H_
#include "cInParam.h"
#include "cHmm.h"
#include "logprob.h"

class cLogBaumWelch
{       private :
                        uint    mvNSample               ;
                        uint*   mvT                             ;
        public :                                                ;
                        cOTMatrix*      mLogAlpha       ;
                        cOTMatrix*      mLogBeta        ;
                        cOTVector*      mLogRho         ;
                        cOTMatrix*      mLogGamma       ;
                        cOTMatrix**     mLogXsi         ;
                        cOTMatrix*      mSumLogXsi      ;
                        cOTVector       mLogVrais       ;
        public :
                cLogBaumWelch(uint theNSample, uint* theT, uint theNClass) ;
                cLogBaumWelch(const cInParam &theInParam) ;
                void LogForwardBackward(cOTMatrix* theCondProba, cHmm& theHMM) ;
                uint GetSampleSize(uint theN){ return mvT[theN] ;}
                virtual ~cLogBaumWelch() ;
} ;


#endif // _CLOGBAUMWELCH_H_
