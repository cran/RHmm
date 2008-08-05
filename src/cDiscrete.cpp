/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
 *                                                         *
 ***********************************************************/
#include "cdiscrete.h"

cDiscrete::cDiscrete(uint theDimObs, uint theNProba) 
{	if ( (theDimObs > 0) && (theNProba > 0) )
	{	mvNClass = theDimObs ;
		mProba = new cOTVector[mvNClass] ;
		for (register uint i = 0 ; i < mvNClass ; i++)
			mProba[i].ReAlloc(theNProba) ;
	}
	else
	{	mvNClass = 0 ; 
		mProba = NULL ;
	}
}

cDiscrete::~cDiscrete()
{	if ( mvNClass > 0)
	{	for (register uint i = 0 ; i < mvNClass ; i++)
			mProba[i].Delete() ;
		delete mProba ;
		mProba = NULL ;
	}
}
uint cDiscrete::GetNProba(void)
{
	if (mvNClass > 0)
		return mProba[0].mSize ;
	else
		return 0 ;
}

void cDiscrete::Print()
{
	for (register uint i = 0 ; i < mvNClass ; i++)
	{	Rprintf("State %d :\t", i) ;
		for (register uint j = 0 ; j < GetNProba() ; j++)
			Rprintf("P[%d]=%lf\t", j, mProba[i][j]) ;
		Rprintf("\n") ;
	}
}

void cDiscrete::ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba)
{
register uint	i,
				n,
				t	;

	for (n = 0 ; n < theNSample ; n++)
		for (i = 0 ; i < mvNClass ; i++)
		{	for (t = 0 ; t < theY[n].mSize ; t++)
				theCondProba[n][i][t] = mProba[i][(uint)theY[n][t]] ;
		}
}

void cDiscrete::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba)
{
register uint	i	;
uint myNProba = GetNProba() ;
	for (i = 0 ; i < mvNClass ; i++)
	{	double myDenominateur = 0.0 ;
		register uint	n,
						t	;
		for (n = 0 ; n < theInParam.mNSample ; n++)
			for (t = 0 ; t < theInParam.mY[n].mSize  ; t++)
				myDenominateur += theBaumWelch.mGamma[n][i][t] ;

		for (register uint k = 0 ; k < myNProba ; k++)
		{	mProba[i][k] = 0.0 ;
			for (n = 0 ; n < theInParam.mNSample ; n++)
				for ( t = 0 ; t < theInParam.mY[n].mSize ; t++)
					mProba[i][k] += theBaumWelch.mGamma[n][i][t]*(theInParam.mY[n][t]==k) ;
			mProba[i][k] /= myDenominateur ;
		}
	}
}

void cDiscrete::InitParameters(cBaumWelchInParam& theInParam)
{
#ifdef _RDLL_
	GetRNGstate();
#endif //_RDLL_

register uint	i			;
uint myNProba = GetNProba() ;
	for (i = 0 ; i < mvNClass ; i++)
	{	register uint	j			;
		double			mySum = 0.0 ;
		for(j = 0 ; j < myNProba ; j++)
		{	mProba[i][j] =  unif_rand() ;
			mySum += mProba[i][j] ;
		}
		mProba[i] /= mySum ;
	}

#ifdef _RDLL_
	PutRNGstate() ;
#endif //_RDLL_
}
void cDiscrete::CopyDistr(cDistribution* theSrc)
{
	cDiscrete *mySrc = (cDiscrete *)theSrc ;

	mvNClass = mySrc->mvNClass ;
	for (register uint i=0 ; i < mvNClass ; i++)
		mProba[i] = mySrc->mProba[i] ;
}


void cDiscrete::GetParam(uint theDeb, cOTVector& theParam)
{
uint myNProba = this->GetNProba();
register uint k = theDeb ;
	for (register uint n = 0 ; n < mvNClass ; n++)
		for (register uint p = 0 ; p < myNProba - 1 ; p++)
			theParam[k++] = mProba[n][p] ;
}
void cDiscrete::SetParam(uint theDeb, cOTVector& theParam)
{
uint myNProba = GetNProba() ;
register uint k = theDeb ;
	for (register uint n = 0 ; n < mvNClass ; n++)
	{	mProba[n][myNProba-1] = 1.0L ;
		for (register uint p = 0 ; p < myNProba - 1 ; p++)
		{	mProba[n][p]  = theParam[k++] ;
			mProba[n][myNProba-1] -= mProba[n][p] ;
		}
	}
}
