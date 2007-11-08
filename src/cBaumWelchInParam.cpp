/***********************************************************
 * RHmm version 0.9.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/11/08                                        *
 *                                                         *
 ***********************************************************/
#include "cBaumWelchInParam.h"


cBaumWelchInParam::cBaumWelchInParam(uint theNSample,uint theDimObs, cOTVector *theY, distrDefinitionEnum theDistrType, uint theNClass, uint theNMixt, uint theNProba, bool theNoHmm):cInParam(theNSample, theDimObs, theY, theDistrType, theNClass, theNMixt, theNProba)
{	
	SetDefault() ;
	mNoHmm = theNoHmm ;
}
cBaumWelchInParam::~cBaumWelchInParam()
{
}
void cBaumWelchInParam::SetDefault(void)
{
	mInitType = eRandom ;
	mNMaxIter = 100 ;
	mTol = 1e-6 ;
	mNInitIter = 5 ;
	mNMaxIterInit = 10 ;
	mVerbose = 0 ;
}

cBaumWelchInParam &cBaumWelchInParam::operator =(const cBaumWelchInParam &theSrc)
{	mInitType = theSrc.mInitType ;
	mNMaxIter = theSrc.mNMaxIter ;
	mTol = theSrc.mTol ;
	mNInitIter = theSrc.mNMaxIterInit ;
	mNMaxIterInit = theSrc.mNMaxIterInit ;
	mVerbose = theSrc.mVerbose ;
	mNoHmm = theSrc.mNoHmm ;

	mDistrType = theSrc.mDistrType ;		
	mNClass = theSrc.mNClass ;
	if (mNSample > 0)
	{	for (register uint i = 0 ; i < mNSample ; i++)
			mY[i].Delete() ;
		delete mY ;
	}
	mY = new cOTVector[theSrc.mNSample] ;
	mNSample = theSrc.mNSample ;
	mDimObs = theSrc.mDimObs ;
	mNProba = theSrc.mNProba ;
	mNMixt = theSrc.mNMixt ;
	for (register uint i = 0 ; i < mNSample ; i++)
		mY[i] = theSrc.mY[i] ;
	return *this ;
}

void cBaumWelchInParam::Print(void)
{
	Rprintf("NbSample = %d\n", mNSample) ;
	for (register uint n = 0 ; n < mNSample ; n++)
		Rprintf("mT[%d]=%d\n", n, (mY[n].mSize)/mDimObs) ;
}
