/**************************************************************
 *** RHmm version 1.2.0                                      
 ***                                                         
 *** File: cHmm.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2008/11/29                                        
 ***                                                         
 **************************************************************/

#include "cHmm.h"
#include "AllDistributions.h"


cHmm::cHmm(distrDefinitionEnum theDistrType, uint theNClass, uint theDimObs, uint theNMixture, uint theNProba)
{	MESS_CREAT("cHmm")
	mDistrType = theDistrType ;
	mInitProba.ReAlloc(theNClass) ;
	mTransMat.ReAlloc(theNClass, theNClass) ;

	switch(mDistrType)
	{	case eNormalDistr :
			mDistrParam = new cUnivariateNormal(theNClass) ;
		break ;
		case eMultiNormalDistr :
			mDistrParam = new cMultivariateNormal(theNClass, theDimObs) ;
		break ;
		case eMixtUniNormalDistr :
			mDistrParam = new cMixtUnivariateNormal(theNClass, theNMixture) ;
		break ;
		case eMixtMultiNormalDistr :
			mDistrParam = new cMixtMultivariateNormal(theNClass, theNMixture, theDimObs) ;
		break ;
		case eDiscreteDistr :
			mDistrParam = new cDiscrete(theNClass, theNProba) ;
		break ;
		case eUnknownDistr :
			mDistrParam = (cDistribution *)NULL ;
		break ;
	}
}
cHmm::cHmm(const cInParam &theInParam)
{	MESS_CREAT("cHmm")
	mInitProba.ReAlloc(theInParam.mNClass);
	mTransMat.ReAlloc(theInParam.mNClass, theInParam.mNClass) ;
	mDistrType = theInParam.mDistrType ;
	switch(mDistrType)
	{	case eNormalDistr :
			mDistrParam = new cUnivariateNormal(theInParam.mNClass) ;
		break ;
		case eMultiNormalDistr :
			mDistrParam = new cMultivariateNormal(theInParam.mNClass, theInParam.mDimObs) ;
		break ;
		case eMixtUniNormalDistr :
			mDistrParam = new cMixtUnivariateNormal(theInParam.mNClass, theInParam.mNMixt) ;
		break ;
		case eMixtMultiNormalDistr :
			mDistrParam = new cMixtMultivariateNormal(theInParam.mNClass, theInParam.mNMixt, theInParam.mDimObs) ;
		break ;
		case eDiscreteDistr :
			mDistrParam = new cDiscrete(theInParam.mNClass, theInParam.mNProba) ;
		break ;
		case eUnknownDistr :
			mDistrParam = (cDistribution *)NULL ;
		break ;
	}
}
cHmm::~cHmm()
{	MESS_DESTR("cHmm")
	mInitProba.Delete() ;
	mTransMat.Delete() ;
	delete mDistrParam ;	
/*	if (mDistrParam != NULL)
	{	switch(mDistrType)
		{	case eNormalDistr :
			{	
			cUnivariateNormal* myDistr = dynamic_cast<cUnivariateNormal *>(mDistrParam) ;
				delete myDistr ;
			}
			break ;
			case eMultiNormalDistr :
			{
			cMultivariateNormal* myDistr = dynamic_cast<cMultivariateNormal *>(mDistrParam) ;
				delete myDistr ;
			}
			break ;
			case eMixtUniNormalDistr :
			{
			cMixtUnivariateNormal* myDistr = dynamic_cast<cMixtUnivariateNormal *>(mDistrParam) ;
				delete myDistr ;
			}
			break ;
			case eMixtMultiNormalDistr :
			{	
			cMixtUnivariateNormal* myDistr = dynamic_cast<cMixtUnivariateNormal *>(mDistrParam) ;
				delete myDistr ;
			}
			break ;
			case eDiscreteDistr :
			{
			cDiscrete* myDistr = dynamic_cast<cDiscrete *>(mDistrParam) ;
				delete myDistr ;
			}
			break ;
			case eUnknownDistr :
				mDistrParam = (cDistribution *)NULL ;
			break ;
		}
	}
*/
	mDistrParam = (cDistribution *)NULL ;
}

cHmm & cHmm::operator = (cHmm &theSrc)
{	/*mvQ = theSrc.mvQ ;
	mvN = theSrc.mvN ;
	mvNMixt = theSrc.mvNMixt ;
	mvNProba = theSrc.mvNProba ;
	*/
	mInitProba = theSrc.mInitProba ;
	mTransMat = theSrc.mTransMat ;

	/*for (register uint i = 0 ; i < mvQ ; i++)
	{	mInitProba[i] = theSrc.mInitProba[i] ;
		for (register uint j = 0 ; j < mvQ ; j++)
			mTransMat[i][j] = theSrc.mTransMat[i][j] ;
	}
	*/
	mDistrParam->CopyDistr(theSrc.mDistrParam) ;
	return(*this) ;
}

void cHmm::Print(void)
{
register uint	i,
				j	;

	Rprintf("ProbInit :\n") ;
	for (i = 0 ; i < mInitProba.mSize ; i++)
		Rprintf("\t%f", mInitProba[i]) ;
	Rprintf("\nMatrice de transition : \n") ;
	for (i = 0 ; i < mInitProba.mSize ; i++)
	{	for (j = 0 ; j < mInitProba.mSize ; j++)
			Rprintf("\t%f", mTransMat[i][j]) ;
		Rprintf("\n") ;
	}
	Rprintf("\nParameters:\n") ;
	mDistrParam->Print() ;
}

uint cHmm::GetNParam(void)
{
uint myNClass = mInitProba.mSize ;
	return( -1 + myNClass * (myNClass + mDistrParam->GetNParam()) ) ;
}

void cHmm::SetParam(cOTVector& theParam) 
{
uint myNClass = mInitProba.mSize ;
register uint k = 0 ;
	
	mInitProba[myNClass-1] = 1.0L ;
	for (register uint n = 0 ; n < myNClass - 1 ; n++)
	{	mInitProba[n] = theParam[k++] ;
		mInitProba[myNClass-1] -= mInitProba[n] ;
	}

	for (register uint n = 0 ; n < myNClass ; n++)
	{	mTransMat[n][myNClass-1] = 1.0L ;
		for (register uint p = 0 ; p < myNClass - 1 ; p++)
		{	mTransMat[n][p] = theParam[k++] ;
			mTransMat[n][myNClass-1] -= mTransMat[n][p] ;
		}
	}
	mDistrParam->SetParam(k, theParam) ;
}

void cHmm::GetParam(cOTVector& theParam) 
{
uint myNClass = mInitProba.mSize ; 
register uint k = 0 ;
	for (register uint n = 0 ; n < myNClass - 1 ; n++)
		theParam[k++] = mInitProba[n] ;
	for (register uint n = 0 ; n < myNClass ; n++)
		for (register uint p = 0 ; p < myNClass - 1 ; p++)
			theParam[k++] = mTransMat[n][p] ;
	mDistrParam->GetParam(k, theParam) ;
}

