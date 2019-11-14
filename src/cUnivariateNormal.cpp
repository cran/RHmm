/**************************************************************
 *** RHmm version 1.3.4                                      
 ***                                                         
 *** File: cUnivariateNormal.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2010/11/14                                      
 ***                                                         
 **************************************************************/

#include "cUnivariateNormal.h"

cUnivariateNormal::cUnivariateNormal(uint theDimObs)
{	MESS_CREAT("cUnivariateNormal")
	mMean.ReAlloc(theDimObs) ;
	mVar.ReAlloc(theDimObs) ;
}

cUnivariateNormal::~cUnivariateNormal()
{	MESS_DESTR("cUnivariateNormal")
	mMean.Delete() ;
	mVar.Delete() ;
}

void cUnivariateNormal::ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba)
{
register uint	i,
				n,
				t						;
double			myAux					;

	for (n = 0 ; n < theNSample ; n++)
		for (i = 0 ; i < mMean.mSize ; i++)
		{	
		double mySigma = sqrt(mVar[i]) ;
			for (t = 0 ; t < theY[n].mSize ; t++)
			{	myAux = (theY[n][t] - mMean[i])/mySigma ;
				theCondProba[n][i][t] = 1.0/(SQRT_TWO_PI*mySigma)*exp(-myAux*myAux/2.0) ;
			}
		}
}
void cUnivariateNormal::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba)
{	
	for (register uint i = 0 ; i < theInParam.mNClass ; i++)
	{	register uint	n,
						t	;
	double myDenominateur = 0.0 ;
		for (n = 0 ; n < theInParam.mNSample ; n++)
			for (t = 0 ; t < theInParam.mY[n].mSize  ; t++)
				myDenominateur += theBaumWelch.mGamma[n][i][t] ;
		mVar[i] = mMean[i] = 0.0 ;
		for (n = 0 ; n < theInParam.mNSample ; n++)
			for (t = 0 ; t < theInParam.mY[n].mSize ; t++)
			{	mMean[i] += theBaumWelch.mGamma[n][i][t] * theInParam.mY[n][t] ;
				mVar[i] += theBaumWelch.mGamma[n][i][t] * theInParam.mY[n][t] * theInParam.mY[n][t] ;
			}
		mMean[i] /= myDenominateur ;
		mVar[i] /= myDenominateur ;
		mVar[i] -= mMean[i] * mMean[i] ;
	}
}

void cUnivariateNormal::InitParameters(cBaumWelchInParam &theInParam)
{
#ifndef _RDLL_
	if (theInParam.mInitType == eKMeans)
	{
	uint myT = 0 ;
	register uint	k	;
		for (k = 0 ; k < theInParam.mNSample ; k++)
			myT += theInParam.mY[k].mSize ;

	int *mySeq = new int[myT], 
		*myNbObs = new int[myT]	;

	cOTVector myY(myT)	;
//		alloc_vecteur(myY, myT) ;
	register uint	s = 0,
					t		;
		for (k = 0 ; k < theInParam.mNSample ; k++)
			for (t = 0 ; t < theInParam.mY[k].mSize ; t++)
				myY[s++] = theInParam.mY[k][t] ;

		KMeans(myY, theInParam.mNClass,  mySeq) ;
	cUnivariateNormal myLoi(theInParam.mNClass) ;
		for (k = 0 ; k < theInParam.mNClass ; k++)
			myLoi.mMean[k] = myLoi.mVar[k] = 0.0 ;
		 
		for (t = 0 ; t < myT ; t++)
		{	k = mySeq[t] ;
			myLoi.mMean[k] = ((double)myNbObs[k]*myLoi.mMean[k] + myY[t])/(double)(myNbObs[k]+1) ;
			myLoi.mVar[k] = ((double)myNbObs[k]*myLoi.mVar[k] + myY[t]*myY[t])/(double)(myNbObs[k]+1) ;
			myNbObs[k]++ ;
		}
		for (k = 0 ; k < theInParam.mNClass ; k++)
			myLoi.mVar[k] -= myLoi.mMean[k] * myLoi.mMean[k] ;
		CopyDistr(&myLoi) ;
		delete mySeq ;
		delete myNbObs ;
		myY.Delete() ;
		return ;
	}
#else
	GetRNGstate();
#endif //_RDLL_

double	myMoy = 0, 
		myVar = 0,
		mystdev		;
double	mys = 0.0L		;		
		for (register uint n = 0 ; n < theInParam.mNSample ; n++)
		{	for (register uint t = 0 ; t < theInParam.mY[n].mSize  ; t++)
			{	myMoy = (mys*myMoy + theInParam.mY[n][t])/(mys+1) ;
				myVar = (mys*myVar + theInParam.mY[n][t]*theInParam.mY[n][t])/(mys+1) ;
				mys++ ;
			}
		}
		myVar -= myMoy*myMoy ;
		mystdev = sqrt(myVar) ;
		for (register uint i = 0 ; i < theInParam.mNClass ; i++)
		{	mMean[i] =  -2*mystdev + myMoy + 2*mystdev * unif_rand() ;
			mVar[i] = 0.5*myVar + 3*myVar * unif_rand() ;	 ;
		}

#ifdef _RDLL_
	PutRNGstate() ;
#endif //_RDLL_
}

void cUnivariateNormal::Print()
{
	Rprintf("Parametres\n") ;
	for (register uint i = 0 ; i < mMean.mSize ; i++)
		Rprintf("m[%d]=%lf\ts[%d]=%f\n", i, mMean[i], i, sqrt(mVar[i]));
}

void cUnivariateNormal::CopyDistr(cDistribution* theSrc)
{
cUnivariateNormal* mySrc = static_cast<cUnivariateNormal *>(theSrc) ;
	mMean = mySrc->mMean ;
	mVar = mySrc->mVar ;
}





void cUnivariateNormal::GetParam(uint theDeb, cOTVector& theParam)
{
register uint k = theDeb ;
	for (register uint n = 0 ; n < mMean.mSize ; n++)
	{	theParam[k++] = mMean[n] ;
		theParam[k++] = mVar[n] ;
	}
}
void cUnivariateNormal::SetParam(uint theDeb, cOTVector& theParam)
{
register uint k = theDeb ;
	for (register uint n = 0 ; n < mMean.mSize ; n++)
	{	mMean[n] = theParam[k++] ;
		mVar[n] = theParam[k++] ;
	}
}

