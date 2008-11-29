/**************************************************************
 *** RHmm version 1.2.0                                      
 ***                                                         
 *** File: cLogBaumWelch.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2008/11/29                                        
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
		mLogRho = NULL ;
		return ;
	}
	mvT = new uint[mvNSample] ;
	mLogVrais.ReAlloc(mvNSample) ;
	
	mLogAlpha = new cOTMatrix[mvNSample] ;
	mLogRho = new cOTVector[mvNSample] ;
	for (register uint n = 0 ; n < mvNSample ; n++)
	{	mvT[n] = theT[n] ;
		mLogAlpha[n].ReAlloc(theNClass, mvT[n]) ;
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
		mLogRho = NULL ;
		return ;
	}	
	mvT = new uint[mvNSample] ;
	mLogVrais.ReAlloc(mvNSample) ;
	
	mLogAlpha = new cOTMatrix[mvNSample] ;
	mLogRho = new cOTVector[mvNSample] ;
	for (register uint n = 0 ; n < mvNSample ; n++)
	{	mvT[n] = (theInParam.mY[n].mSize)/theInParam.mDimObs ;
		mLogAlpha[n].ReAlloc(theInParam.mNClass, mvT[n]) ;
		mLogRho[n].ReAlloc(mvT[n]) ;
	}	
}

cLogBaumWelch::~cLogBaumWelch()
{	MESS_DESTR("cLogBaumWelch") 
	if (mvNSample > 0)
	{	for (register uint n = 0 ; n < mvNSample ; n++)
		{	mLogAlpha[n].Delete() ;
			mLogRho[n].Delete() ;
		}
		delete [] mvT ;
		delete [] mLogRho ;
		delete [] mLogAlpha ;
	}
}


void cLogBaumWelch::LogForwardBackward(cOTMatrix* theCondProba, cHmm& theHMM)
{
register uint	i,
				j		;
register int	t		;
double			myLogAlpha ;
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
		for ( i = 0 ; i < myNClass ; i++)
			mLogAlpha[n][i][0] = elnproduct(mLogAlpha[n][i][0], -mLogRho[n][0]) ; // Normalisation
	//forward

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
			for (j = 0 ; j < myNClass ; j++)
				mLogAlpha[n][j][t+1] = elnproduct(mLogAlpha[n][j][t+1], -mLogRho[n][t+1]) ;

			mLogVrais[n] = elnproduct(mLogVrais[n], mLogRho[n][t+1]) ;
		}
	}
}
