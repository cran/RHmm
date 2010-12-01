/**************************************************************
 *** RHmm version 1.4.3                                     
 ***                                                         
 *** File: cBaumWelch.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/01                                     
 ***                                                         
 **************************************************************/

#include "cBaumWelch.h"


cBaumWelch::cBaumWelch(uint theNSample, uint* theT, uint theNClass)
{       MESS_CREAT("cBaumWelch") 
        mvNSample = theNSample ;
        if (mvNSample == 0)
        {       mvT = NULL ;
                mLogVrais.Delete() ;
                mAlpha = NULL ;
                mBeta = NULL ;
                mGamma = NULL ;
                mXsi = NULL ;
                mSumXsi = NULL ;
                mRho = NULL ;
                return ;
        }
        mvT = new uint[mvNSample] ;
        mLogVrais.ReAlloc(mvNSample) ;
        
        mAlpha = new cOTMatrix[mvNSample] ;
        mBeta = new cOTMatrix[mvNSample] ;
        mGamma = new cOTMatrix[mvNSample] ;
        mXsi = new cOTMatrix*[mvNSample] ;
        mSumXsi = new cOTMatrix[mvNSample] ;
        mRho = new cOTVector[mvNSample] ;
        for (register uint n = 0 ; n < mvNSample ; n++)
        {       mvT[n] = theT[n] ;
                mAlpha[n].ReAlloc(theNClass, mvT[n]) ;
                mBeta[n].ReAlloc(theNClass, mvT[n]) ;
                mGamma[n].ReAlloc(theNClass, mvT[n]) ;
                mXsi[n] = new cOTMatrix[mvT[n]] ;
                for (register uint t = 0 ; t < mvT[n] ; t++)
                        mXsi[n][t].ReAlloc(theNClass, theNClass) ;
                mSumXsi[n].ReAlloc(theNClass, theNClass) ;
                mRho[n].ReAlloc(mvT[n]) ;
        }       
}

cBaumWelch::cBaumWelch(const cInParam &theInParam)
{       MESS_CREAT("cBaumWelch") 
        mvNSample = theInParam.mNSample ;
        if (mvNSample == 0)
        {       mvT = NULL ;
                mLogVrais.Delete() ;
                mAlpha = NULL ;
                mBeta = NULL ;
                mGamma = NULL ;
                mXsi = NULL ;
                mRho = NULL ;
                return ;
        }       
        mvT = new uint[mvNSample] ;
        mLogVrais.ReAlloc(mvNSample) ;
        
        mAlpha = new cOTMatrix[mvNSample] ;
        mBeta = new cOTMatrix[mvNSample] ;
        mGamma = new cOTMatrix[mvNSample] ;
        mXsi = new cOTMatrix*[mvNSample] ;
        mSumXsi = new cOTMatrix[mvNSample] ;
        mRho = new cOTVector[mvNSample] ;
        for (register uint n = 0 ; n < mvNSample ; n++)
        {       mvT[n] = (theInParam.mY[n].mSize)/theInParam.mDimObs ;
                mAlpha[n].ReAlloc(theInParam.mNClass, mvT[n]) ;
                mBeta[n].ReAlloc(theInParam.mNClass, mvT[n]) ;
                mGamma[n].ReAlloc(theInParam.mNClass, mvT[n]) ;
                mXsi[n] = new cOTMatrix[mvT[n]] ;
                for (register uint t=0 ; t < mvT[n] ; t++)
                        mXsi[n][t].ReAlloc(theInParam.mNClass, theInParam.mNClass) ;
                mSumXsi[n].ReAlloc(theInParam.mNClass, theInParam.mNClass) ;
                mRho[n].ReAlloc(mvT[n]) ;
        }       
}

cBaumWelch::~cBaumWelch()
{       MESS_DESTR("cBaumWelch") 
        if (mvNSample > 0)
        {       for (register uint n = 0 ; n < mvNSample ; n++)
                {       mAlpha[n].Delete() ;
                        mBeta[n].Delete() ;
                        mGamma[n].Delete() ;
                        for (register uint t = 0 ; t < mvT[n] ; t++)
                                mXsi[n][t].Delete() ;
                        delete [] mXsi[n] ;
                        mSumXsi[n].Delete() ;
                        mRho[n].Delete() ;
                }
                delete [] mvT ;
                delete [] mRho ;
                delete [] mXsi ;
                delete [] mSumXsi ;
                delete [] mGamma ;
                delete [] mBeta ;
                delete [] mAlpha ;
        }
}

void cBaumWelch::ForwardBackward(cOTMatrix* theCondProba, cHmm& theHMM)
{
register uint   i,
                                j               ;
register int    t               ;
double                  myAux,
                                mySum   ;
uint myNClass = theHMM.mInitProba.mSize ;
        
        for (register uint n = 0 ; n < mvNSample ; n++)
        {
        int myT = (int)mvT[n] ;
                mRho[n][0] = 0.0L ;
                for (i = 0 ; i < myNClass ; i++)
                {       mAlpha[n][i][0] = theHMM.mInitProba[i] * theCondProba[n][i][0] ;
                        mRho[n][0] += mAlpha[n][i][0] ; 
                }
                for ( i = 0 ; i < myNClass ; i++)
                        mAlpha[n][i][0] /= mRho[n][0] ; // Normalisation
        //forward
                for (t = 0 ; t < myT-1 ; t++)
                {       mRho[n][t+1] = 0.0 ;
                        for (j = 0 ; j < myNClass ; j++)
                        {       myAux = 0.0 ;
                                for (i = 0 ; i < myNClass ; i++)
                                        myAux += mAlpha[n][i][t] * theHMM.mTransMatVector[t][i][j] ; /* FIXME: Is t correct or should we shift? */
                                mAlpha[n][j][t+1] = myAux * theCondProba[n][j][t+1] ;
                                mRho[n][t+1] += mAlpha[n][j][t+1] ;
                        }
                        for (j = 0 ; j < myNClass ; j++)
                                mAlpha[n][j][t+1] /= mRho[n][t+1] ;
                }

        // backward
                for (i = 0 ; i < myNClass ; i++)
                        mBeta[n][i][myT-1] = 1.0/mRho[n][myT-1] ;

                for (t = myT-2 ; t >= 0 ; t--)
                {       for (i = 0 ; i < myNClass ; i++)
                        {       myAux = 0.0 ;
                                for (j = 0 ; j < myNClass ; j++)
                                        myAux +=  theHMM.mTransMatVector[t][i][j] * theCondProba[n][j][t+1] * mBeta[n][j][t+1] ; /* FIXME: Is t correct or should we shift? */
                                mBeta[n][i][t] = myAux ;
                        }
                        for (i = 0 ; i < myNClass ; i++)
                                mBeta[n][i][t] /= mRho[n][t] ;
                }
                
        // Calcul des Gamma et des Xsi et LogVrais
                mLogVrais[n] = 0 ;
                for (t = 0 ; t < myT ; t++)
                {       mySum = 0.0 ;
                        for (i = 0 ; i < myNClass ; i++)
                        {       mGamma[n][i][t] = mAlpha[n][i][t] * mBeta[n][i][t] ;
                                mySum += mGamma[n][i][t] ;
                        }
                        for (i = 0 ; i < myNClass ; i++)
                                mGamma[n][i][t] /= mySum ;
                        
                        mLogVrais[n] += log(mRho[n][t]) ;
                }
        // Calcul des Xsi
                for (i = 0 ; i < myNClass ; i++)
                        for (j = 0 ; j < myNClass ; j++)
                        {       mSumXsi[n][i][j] = 0.0 ;
                                for (t = 0 ; t < myT - 1 ; t++)
                                {       mXsi[n][t][i][j] = mAlpha[n][i][t] * theHMM.mTransMatVector[t][i][j] * theCondProba[n][j][t+1] * mBeta[n][j][t+1] ;
                                        mSumXsi[n][i][j] += mXsi[n][t][i][j] ;
                                }
                        }
        }
}


