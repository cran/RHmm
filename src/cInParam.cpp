/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
 *                                                         *
 ***********************************************************/
#include "cinparam.h"

cInParam::cInParam(uint theNSample, uint theDimObs, cOTVector* theY, distrDefinitionEnum theDistrType, uint theNClass, uint theNMixt, uint theNProba)
{	mDistrType = theDistrType ;
	mNClass = theNClass ;
	mNMixt = theNMixt ;
	mNProba = theNProba ;
	mNSample = theNSample ;
	mDimObs = theDimObs ;
	if (mNSample > 0)
	{	mY = new cOTVector[mNSample] ;
		for (register uint i = 0 ; i < mNSample ; i++)
			mY[i] = theY[i] ;
	}
	else
		mY = (cOTVector *)NULL ;
}
cInParam::~cInParam()
{	if (mNSample != 0)
	{	for (register uint i = 0 ; i < mNSample ; i++)
			mY[i].Delete() ;
		delete [] mY ;
		mNSample = 0 ;
	}
}

cInParam &cInParam::operator =(const cInParam &theSrc)
{	
	mDistrType = theSrc.mDistrType ;		
	mNClass = theSrc.mNClass ;
	if (mNSample > 0)
	{	for (register uint i = 0 ; i < mNSample ; i++)
			mY[i].Delete() ;
		delete mY ;
	}
	mNSample = theSrc.mNSample ;
	mY = new cOTVector[mNSample] ;
	
	mDimObs = theSrc.mDimObs ;
	mNProba = theSrc.mNProba ;
	mNMixt = theSrc.mNMixt ;
	
	for (register uint i = 0 ; i < mNSample ; i++)
		mY[i] = theSrc.mY[i] ;
	
	return(*this) ;
}

void cInParam::Print(void)
{
	Rprintf("NbSample = %d\n", mNSample) ;
	for (register uint n = 0 ; n < mNSample ; n++)
		Rprintf("mT[%d]=%d\n", n, mY[n].mSize) ;
}
