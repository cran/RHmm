/**************************************************************
 *** RHmm version 1.3.0                                      
 ***                                                         
 *** File: cLogBaumWelch.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2009/06/03                                      
 ***                                                         
 **************************************************************/

#include "cLogBaumWelch.h"


cLogBaumWelch::cLogBaumWelch(uint theNSample, uint* theT, uint theNClass)
{	MESS_CREAT("cLogBaumWelch") 
	mvNSample = theNSample ;
	if (mvNSample == 0)
	{	mvT = NULL ;
		mLogVrais.Delete() ;
		mLogAlpha = NULL ;
		mLogBeta = NULL ;
		mLogGamma = NULL ;
		mLogXsi = NULL ;
		mSumLogXsi = NULL ;
		mLogRho = NULL ;
		return ;
	}
	mvT = new uint[mvNSample] ;
	mLogVrais.ReAlloc(mvNSample) ;
	
	mLogAlpha = new cOTMatrix[mvNSample] ;
	mLogBeta = new cOTMatrix[mvNSample] ;
	mLogGamma = new cOTMatrix[mvNSample] ;
	mLogXsi = new cOTMatrix*[mvNSample] ;
	mSumLogXsi = new cOTMatrix[mvNSample] ;
	mLogRho = new cOTVector[mvNSample] ;
	for (register uint n = 0 ; n < mvNSample ; n++)
	{	mvT[n] = theT[n] ;
		mLogAlpha[n].ReAlloc(theNClass, mvT[n]) ;
		mLogBeta[n].ReAlloc(theNClass, mvT[n]) ;
		mLogGamma[n].ReAlloc(theNClass, mvT[n]) ;
		mLogXsi[n] = new cOTMatrix[mvT[n]] ;
		for (register uint t = 0 ; t < mvT[n] ; t++)
			mLogXsi[n][t].ReAlloc(theNClass, theNClass) ;
		mSumLogXsi[n].ReAlloc(theNClass, theNClass) ;
		mLogRho[n].ReAlloc(mvT[n]) ;
	}	
}

cLogBaumWelch::cLogBaumWelch(const cInParam &theInParam)
{	MESS_CREAT("cLogBaumWelch") 
	mvNSample = theInParam.mNSample ;
	if (mvNSample == 0)
	{	mvT = NULL ;
		mLogVrais.Delete() ;
		mLogAlpha = NULL ;
		mLogBeta = NULL ;
		mLogGamma = NULL ;
		mLogXsi = NULL ;
		mLogRho = NULL ;
		return ;
	}	
	mvT = new uint[mvNSample] ;
	mLogVrais.ReAlloc(mvNSample) ;
	
	mLogAlpha = new cOTMatrix[mvNSample] ;
	mLogBeta = new cOTMatrix[mvNSample] ;
	mLogGamma = new cOTMatrix[mvNSample] ;
	mLogXsi = new cOTMatrix*[mvNSample] ;
	mSumLogXsi = new cOTMatrix[mvNSample] ;
	mLogRho = new cOTVector[mvNSample] ;
	for (register uint n = 0 ; n < mvNSample ; n++)
	{	mvT[n] = (theInParam.mY[n].mSize)/theInParam.mDimObs ;
		mLogAlpha[n].ReAlloc(theInParam.mNClass, mvT[n]) ;
		mLogBeta[n].ReAlloc(theInParam.mNClass, mvT[n]) ;
		mLogGamma[n].ReAlloc(theInParam.mNClass, mvT[n]) ;
		mLogXsi[n] = new cOTMatrix[mvT[n]] ;
		for (register uint t=0 ; t < mvT[n] ; t++)
			mLogXsi[n][t].ReAlloc(theInParam.mNClass, theInParam.mNClass) ;
		mSumLogXsi[n].ReAlloc(theInParam.mNClass, theInParam.mNClass) ;
		mLogRho[n].ReAlloc(mvT[n]) ;
	}	
}

cLogBaumWelch::~cLogBaumWelch()
{	MESS_DESTR("cLogBaumWelch") 
	if (mvNSample > 0)
	{	for (register uint n = 0 ; n < mvNSample ; n++)
		{	mLogAlpha[n].Delete() ;
			mLogBeta[n].Delete() ;
			mLogGamma[n].Delete() ;
			for (register uint t = 0 ; t < mvT[n] ; t++)
				mLogXsi[n][t].Delete() ;
			delete [] mLogXsi[n] ;
			mSumLogXsi[n].Delete() ;
			mLogRho[n].Delete() ;
		}
		delete [] mvT ;
		delete [] mLogRho ;
		delete [] mLogXsi ;
		delete [] mSumLogXsi ;
		delete [] mLogGamma ;
		delete [] mLogBeta ;
		delete [] mLogAlpha ;
	}
}

void cLogBaumWelch::LogForwardBackward(cOTMatrix* theCondProba, cHmm& theHMM)
{
register uint	i,
				j		;
register int	t		;
double			myLogAlpha,
				myAux,
				mySum	;
uint myNClass = theHMM.mInitProba.mSize ;
cOTMatrix myLogTransMat = cOTMatrix(myNClass, myNClass) ;
cOTVector myLogInitProba = cOTVector(myNClass) ;
	
	for ( i = 0 ; i < myNClass ; i++)
	{	myLogInitProba[i] = eln(theHMM.mInitProba[i]) ;
		for (j = 0 ; j < myNClass ; j++)
			myLogTransMat[i][j] = eln(theHMM.mTransMat[i][j]) ;
	}

	for (register uint n = 0 ; n < mvNSample ; n++)
	{
	int myT = (int)mvT[n] ;
		mLogRho[n][0] = LOGZERO ;
		for (i = 0 ; i < myNClass ; i++)
		{	mLogAlpha[n][i][0] = elnproduct(myLogInitProba[i], eln(theCondProba[n][i][0])) ;
			mLogRho[n][0] = elnsum(mLogRho[n][0], mLogAlpha[n][i][0]) ;	
		}

		mLogVrais[n] = mLogRho[n][0] ;
		//forward
		for (t = 0 ; t < myT-1 ; t++)
		{	mLogRho[n][t+1] = LOGZERO ;
			for (j = 0 ; j < myNClass ; j++)
			{	myLogAlpha = LOGZERO ;
				for (i = 0 ; i < myNClass ; i++)
					myLogAlpha = elnsum(myLogAlpha, elnproduct(mLogAlpha[n][i][t], myLogTransMat[i][j])) ;
				
				mLogAlpha[n][j][t+1] = elnproduct(myLogAlpha, eln(theCondProba[n][j][t+1])) ;
				mLogRho[n][t+1] = elnsum(mLogRho[n][t+1], mLogAlpha[n][j][t+1]) ;
			}
		}

		// backward
		for (i = 0 ; i < myNClass ; i++)
			mLogBeta[n][i][myT-1] = 0.0L ;

		for (t = myT-2 ; t >= 0 ; t--)
		{	for (i = 0 ; i < myNClass ; i++)
			{	myAux = LOGZERO ;
				for (j = 0 ; j < myNClass ; j++)
					myAux =  elnsum(myAux, elnproduct(theHMM.mTransMat[i][j], elnproduct(theCondProba[n][j][t+1], mLogBeta[n][j][t+1]))) ;
				mLogBeta[n][i][t] = myAux ;
			}
		}
		
		// Calcul des Gamma et LogVrais
		mLogVrais[n] = LOGZERO ;
		for (t = 0 ; t < myT ; t++)
		{	mySum = LOGZERO ;
			for (i = 0 ; i < myNClass ; i++)
			{	mLogGamma[n][i][t] = elnproduct(mLogAlpha[n][i][t], mLogBeta[n][i][t]) ;
				mySum = elnsum(mySum, mLogGamma[n][i][t]) ;
			}
			for (i = 0 ; i < myNClass ; i++)
				mLogGamma[n][i][t] = elnproduct(mLogGamma[n][i][t], -mySum) ;
			
			mLogVrais[n] = elnsum(mLogVrais[n], mLogRho[n][t]) ;
		}
	
	// Calcul des Xsi
		for (i = 0 ; i < myNClass ; i++)
			for (j = 0 ; j < myNClass ; j++)
			{	mSumLogXsi[n][i][j] = LOGZERO ;
				for (t = 0 ; t < myT - 1 ; t++)
				{	mLogXsi[n][t][i][j] = elnproduct(mLogAlpha[n][i][t], elnproduct(theHMM.mTransMat[i][j], elnproduct(theCondProba[n][j][t+1], mLogBeta[n][j][t+1]))) ;
					mSumLogXsi[n][i][j] = elnsum(mSumLogXsi[n][i][j],mLogXsi[n][t][i][j]) ;
				}
			}

	// fin
		for (i = 0 ; i < myNClass ; i++)
		{	for (t = 0 ; t < myT ; t++)
			{	mLogAlpha[n][i][t] = eexp(mLogAlpha[n][i][t]) ;
				mLogBeta[n][i][t] = eexp(mLogBeta[n][i][t]) ;
				mLogGamma[n][i][t] = eexp(mLogGamma[n][i][t]) ;
			}
			for (j = 0 ; j < myNClass ; j++)
			{	mSumLogXsi[n][i][j] = eexp(mSumLogXsi[n][i][j]) ;
				for (t = 0 ; t < myT - 1 ; t++)
					mLogXsi[n][t][i][j] = eexp(mLogXsi[n][t][i][j]) ;
			}
		}
		for (t = 0 ; t < myT ; t++)
			mLogRho[n][t] = eexp(mLogRho[n][t]) ;

	}
}
