/**************************************************************
 *** RHmm version 1.3.0                                      
 ***                                                         
 *** File: RHmm.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2009/04/27                                       
 ***                                                         
 **************************************************************/

#ifdef _RDLL_

#include "OTMathUtil.h"
#include "cRUtils.h"
#include "Hmm.h"
#include "cInParam.h"
#include "cBaumWelchInParam.h"
#include "cBaumWelch.h"
#include "cLogBaumWelch.h"
#include "cLogBaumWelch.h"
#include "AllDistributions.h"
#include "cHmm.h"
#include "cHmmFit.h"
#include "cViterbi.h"
#include "RHmm.h"


BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RBaumWelch	(	SEXP theParamHMM, 
									SEXP theYt, 
									SEXP theParamBW
								)
{
initEnum		myTypeInit			;
distrDefinitionEnum	myDistrType		;
uint				myNbClasses,
					myDimObs,
					myNbMixt,
					myNbProba		;
cRUtil				myRUtil			;
// Récupération des paramètres d'entrée
	myRUtil.GetValSexp(theParamHMM, eNClasses, myNbClasses) ;
	myRUtil.GetValSexp(theParamHMM, eObsDim, myDimObs) ;
	myRUtil.GetValSexp(theParamHMM, eNMixt, myNbMixt) ;
	myRUtil.GetValSexp(theParamHMM, eNProba, myNbProba) ;

char myStr[255] ;
char *myString = (char *)myStr ;
	myRUtil.GetValSexp(theParamHMM, eDistrType, myString) ;
	myDistrType = eUnknownDistr ;
	if (strcmp(myString, "NORMAL") == 0)
	{		if (myDimObs == 1)
				myDistrType = eNormalDistr ;
			else
				myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myString, "DISCRETE") == 0)
				myDistrType = eDiscreteDistr ;
		else
		{	if (strcmp(myString, "MIXTURE") == 0)
				if (myDimObs==1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
		}
	}

	myRUtil.GetValSexp(theParamBW, eInitType, (char *)myString) ;
	if (strcmp(myString, "RANDOM") == 0)
		myTypeInit = eRandom ;
	else
	{	if (strcmp(myString, "KMEANS") == 0)
			myTypeInit = eKMeans ;
		else
			if (strcmp(myString, "USER") == 0)
				myTypeInit = eUser ;
	}

uint myNbIterMax ;
	myRUtil.GetValSexp(theParamBW, eNMaxIter, myNbIterMax) ;
double myTol ;
	myRUtil.GetValSexp(theParamBW, eTol, myTol) ;
uint myVerbose ;
	myRUtil.GetValSexp(theParamBW, eVerbose, myVerbose) ;
uint myNbIterInit ;
	myRUtil.GetValSexp(theParamBW, eNInitIter, myNbIterInit) ;
uint myNbIterMaxInit ;
	myRUtil.GetValSexp(theParamBW, eNMaxIterinit, myNbIterMaxInit) ;

uint myNbSample = length(theYt) ;

cOTVector* myY = new cOTVector[myNbSample] ;

	for (register uint n = 0 ; n < myNbSample ; n++)
	{	
	SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myY[n].ReAlloc(length(myAux), REAL(myAux)) ;
	}

cBaumWelchInParam myParamEntree = cBaumWelchInParam(myNbSample, myDimObs, myY, myDistrType, myNbClasses, myNbMixt, myNbProba) ;

	myParamEntree.mNMaxIter = myNbIterMax ;
	myParamEntree.mTol = myTol ;
	myParamEntree.mVerbose = myVerbose ;
	myParamEntree.mNInitIter = myNbIterInit ;
	myParamEntree.mNMaxIterInit = myNbIterMaxInit ;
	myParamEntree.mInitType = myTypeInit ;	

cHmmFit myParamSortie = cHmmFit(myParamEntree) ;

	if (myTypeInit == eUser)
	{	cHmm myHMM = cHmm(myParamEntree) ;
	SEXP myAux1 ;
		myRUtil.GetValSexp(theParamBW, eInitPoint, myAux1) ;
		myRUtil.GetVectSexp(myAux1, 0, myHMM.mInitProba) ;
		myRUtil.GetMatSexp(myAux1, 1, myHMM.mTransMat) ;
		SEXP myAux ;
		myRUtil.GetValSexp(myAux1, 2, myAux) ; // $distribution
		switch (myDistrType)
		{	case eNormalDistr :
			{	cUnivariateNormal* myParam = dynamic_cast<cUnivariateNormal *>(myHMM.mDistrParam) ;			
				myRUtil.GetVectSexp(myAux, 3, myParam->mMean) ;
				myRUtil.GetVectSexp(myAux, 4, myParam->mVar) ;
			}
			break ;
			case eMultiNormalDistr :
			{	cMultivariateNormal *myParam = dynamic_cast<cMultivariateNormal *>(myHMM.mDistrParam) ;
				myRUtil.GetListVectSexp(myAux, 3, myNbClasses, myParam->mMean) ; 
				myRUtil.GetListMatSexp(myAux, 4, myNbClasses, myParam->mCov) ;
			}
			break ;
			case eMixtUniNormalDistr :
			{	cMixtUnivariateNormal *myParam = dynamic_cast<cMixtUnivariateNormal *>(myHMM.mDistrParam) ;
				myRUtil.GetListVectSexp(myAux, 4, myNbClasses, myParam->mMean) ;
				myRUtil.GetListVectSexp(myAux, 5, myNbClasses, myParam->mVar) ;
				myRUtil.GetListVectSexp(myAux, 6, myNbClasses, myParam->mp) ;
			}
			break ;
			case eMixtMultiNormalDistr :
			{	cMixtMultivariateNormal *myParam = dynamic_cast<cMixtMultivariateNormal *>(myHMM.mDistrParam) ;
				myRUtil.GetListListVectSexp(myAux, 4, myNbClasses, myNbMixt, myParam->mMean) ;
				myRUtil.GetListListMatSexp(myAux, 5, myNbClasses, myNbMixt, myParam->mCov) ;
				myRUtil.GetListVectSexp(myAux, 6, myNbClasses, myParam->mp) ;
			}
			break ;
			case eDiscreteDistr :
			{	cDiscrete *myParam = dynamic_cast<cDiscrete *>(myHMM.mDistrParam) ;
				myRUtil.GetListVectSexp(myAux, 3, myNbClasses, myParam->mProba) ;
			}
			break ;
		}
		myParamSortie.CopyHmm(myHMM) ;
	}
	else
		myParamSortie.BaumWelchAlgoInit(myParamEntree) ;	
	myParamSortie.BaumWelchAlgo(myParamEntree) ;
	
	for (register uint n = 0 ; n < myNbSample ; n++)
		myY[n].Delete() ;
	delete[] myY ;

/*
#ifdef _DEBOGAGE_
SEXP myRes1 ;
SEXP myAux1[6] ;
uint ourNProtect = 0 ;
	MESS("JE SUIS ICI 1\n") 
	myParamSortie.Print() ;
	SETVECTSEXP(myParamSortie.mInitProba, myAux1[0]) ;	
	SETMATSEXP(myParamSortie.mTransMat, myAux1[1]) ;
cDiscrete *myParam = dynamic_cast<cDiscrete *>(myParamSortie.mDistrParam) ;
	if (myParam == NULL)
		cOTError("Probleme de casting") ;
	else
		MESS("cDiscrete Fait\n") ;
	SETLISTVECTSEXP(myParam->mProba, myNbClasses, myAux1[2]) ;
	PROTECT(myAux1[5] = allocVector(REALSXP, 4)) ;
	REAL(myAux1[5])[0] = myParamSortie.mLLH ;
	REAL(myAux1[5])[1] = myParamSortie.mBic ;
	REAL(myAux1[5])[2] = (double)( myParamSortie.mNIter) ;
	REAL(myAux1[5])[3] =  myParamSortie.mTol ;
	PROTECT(myRes1 = allocVector(VECSXP, 4)) ;
	for (register uint i = 0 ; i < 3 ; i++)
		SET_VECTOR_ELT(myRes1, i, myAux1[i]) ;
	SET_VECTOR_ELT(myRes1, 3, myAux1[5]) ;
	MESS("JE SUIS ICI 2\n") 
	UNPROTECT(2 + ourNProtect) ;
	MESS("JE SUIS ICI 3\n") 
	return(myRes1);
#endif //DEBOGAGE_
*/
SEXP myRes,
	 myAux[6]	;

	myRUtil.SetVectSexp(myParamSortie.mInitProba, myAux[0]) ;
	myRUtil.SetMatSexp(myParamSortie.mTransMat, myAux[1]) ;
	switch (myDistrType)
	{	case eNormalDistr :
		{	cUnivariateNormal* myParam =  dynamic_cast<cUnivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetVectSexp(myParam->mMean, myAux[2]) ;
			myRUtil.SetVectSexp(myParam->mVar, myAux[3]) ;
		}
		break ;
		case eMultiNormalDistr :
		{	cMultivariateNormal *myParam = dynamic_cast<cMultivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListVectSexp(myParam->mMean, myNbClasses, myAux[2]) ;
			myRUtil.SetListMatSexp(myParam->mCov, myNbClasses, myAux[3]) ;
		}
		break ;
		case eDiscreteDistr :
		{	cDiscrete *myParam = dynamic_cast<cDiscrete *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListVectSexp(myParam->mProba, myNbClasses, myAux[2]) ;
		}
		break ;
		case eMixtUniNormalDistr :
		{	cMixtUnivariateNormal *myParam = dynamic_cast<cMixtUnivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListVectSexp(myParam->mMean, myNbClasses, myAux[2]) ;
			myRUtil.SetListVectSexp(myParam->mVar, myNbClasses, myAux[3]) ;
			myRUtil.SetListVectSexp(myParam->mp, myNbClasses,myAux[4]) ;
		}
		break ;
		case eMixtMultiNormalDistr :
		{	cMixtMultivariateNormal *myParam = dynamic_cast<cMixtMultivariateNormal *>(myParamSortie.mDistrParam) ;
			myRUtil.SetListListVectSexp(myParam->mMean, myNbClasses, myNbMixt, myAux[2]) ;
			myRUtil.SetListListMatSexp(myParam->mCov, myNbClasses, myNbMixt, myAux[3]) ;
			myRUtil.SetListVectSexp(myParam->mp, myNbClasses, myAux[4]) ;
		}
		break ;
		default :
		break ;
	}
	PROTECT(myAux[5] = allocVector(REALSXP, 4)) ;
	REAL(myAux[5])[0] = myParamSortie.mLLH ;
	REAL(myAux[5])[1] = myParamSortie.mBic ;
	REAL(myAux[5])[2] = (double)( myParamSortie.mNIter) ;
	REAL(myAux[5])[3] =  myParamSortie.mTol ;

	switch (myDistrType)
	{	case eNormalDistr :
		case eMultiNormalDistr :
			PROTECT(myRes = allocVector(VECSXP, 5)) ;
		break ;
		case eMixtUniNormalDistr :
		case eMixtMultiNormalDistr :
			PROTECT(myRes = allocVector(VECSXP, 6)) ;
		break ;
		case eDiscreteDistr :
			PROTECT(myRes = allocVector(VECSXP, 4)) ;
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

	for (register uint i = 0 ; i < 3 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	switch (myDistrType)
	{	case eNormalDistr :
		case eMultiNormalDistr :
			SET_VECTOR_ELT(myRes, 3, myAux[3]) ;
			SET_VECTOR_ELT(myRes, 4, myAux[5]) ;
		break ;
		case eMixtUniNormalDistr :
		case eMixtMultiNormalDistr :
			SET_VECTOR_ELT(myRes, 3, myAux[3]) ;
			SET_VECTOR_ELT(myRes, 4, myAux[4]) ;
			SET_VECTOR_ELT(myRes, 5, myAux[5]) ;
		break ;
		case eDiscreteDistr :
			SET_VECTOR_ELT(myRes, 3, myAux[5]) ;
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

	UNPROTECT(2) ;
	myRUtil.EndProtect() ;
	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RViterbi	(	SEXP theHMM, 
									SEXP theYt
								)
{
distrDefinitionEnum		myDistrType		;
uint					myDimObs=1,
						myNbClasses,
						myNbProba=0,
						myNbMixt=0		;
cRUtil					myRUtil			;

SEXP myDistSEXP ;
	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba	
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNbClasses) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{	myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
		if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myStr, "DISCRETE") == 0)
		{	myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNbProba) ;
		}
		else
		{	if (strcmp(myStr, "MIXTURE") == 0)
			{	myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
				myRUtil.GetValSexp(myDistSEXP, 2, myNbMixt) ;
			}
		}
	}
uint	myNbSample = length(theYt) ;	
uint*	myT = new uint[myNbSample]	;
//double	**myY	;
cOTVector* myY = new cOTVector[myNbSample] ;
for (register uint n = 0 ; n < myNbSample ; n++)
	{	SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm	myHMM = cHmm(myDistrType, myNbClasses, myDimObs, myNbMixt, myNbProba) ;

	myRUtil.GetVectSexp(theHMM, fInitProba, myHMM.mInitProba) ;
	myRUtil.GetMatSexp(theHMM, fTransMat, myHMM.mTransMat) ;
	
	switch (myDistrType)
	{	case eNormalDistr :
		{	cUnivariateNormal *myLoi = (cUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{	cMultivariateNormal *myLoi = (cMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNbClasses, myLoi->mCov) ;
		}
		break ;
		case  eMixtUniNormalDistr :
		{	cMixtUnivariateNormal *myParam = (cMixtUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNbClasses, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNbClasses, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;
		case  eMixtMultiNormalDistr :
		{	cMixtMultivariateNormal *myParam = (cMixtMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNbClasses, myNbMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNbClasses, myNbMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{	cDiscrete *myParam = (cDiscrete *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myParam->mProba) ;
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}
				
cInParam myParamEntree(myNbSample, myDimObs, myY) ;
	myParamEntree.mDimObs = myDimObs ;
	myParamEntree.mNMixt = myNbMixt ;
	myParamEntree.mNProba = myNbProba ;
	myParamEntree.mNClass = myNbClasses ;
	myParamEntree.mDistrType = myDistrType ;
cViterbi myViterbi = cViterbi(myParamEntree) ;
	myViterbi.ViterbiPath(myParamEntree, myHMM) ;

SEXP myAux[2] ;
	myRUtil.SetListVectSexp(myViterbi.mSeq, myNbSample, myT, myAux[0]) ;
	myRUtil.SetListValSexp(myViterbi.mLogProb, myAux[1]) ;

	SEXP myRes ;
	PROTECT(myRes = allocVector(VECSXP, 2)) ;
	for (register uint i = 0 ; i < 2 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	myRUtil.EndProtect() ;
	UNPROTECT(1) ;
	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP Rforwardbackward	(	SEXP theHMM, 
											SEXP theYt
										)
{
distrDefinitionEnum		myDistrType	;
uint					myDimObs=1,
						myNbClasses,
						myNbProba=0,
						myNbMixt=0	;
cRUtil					myRUtil			;


SEXP myDistSEXP ;

	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba	
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNbClasses) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{	myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
		if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myStr, "DISCRETE") == 0)
		{	myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNbProba) ;
		}
		else
		{	if (strcmp(myStr, "MIXTURE") == 0)
			{	myRUtil.GetValSexp(myDistSEXP, 2, myNbMixt) ;
				myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
			}
		}
	}
uint	myNbSample = length(theYt) ;	
uint*	myT = new uint[myNbSample] ;

cOTVector* myY = new cOTVector[myNbSample] ;


	for (register uint n = 0 ; n < myNbSample ; n++)
	{	SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm myHMM = cHmm(myDistrType, myNbClasses, myDimObs, myNbMixt, myNbProba) ;
	myRUtil.GetVectSexp(theHMM, fInitProba, myHMM.mInitProba) ;
	myRUtil.GetMatSexp(theHMM, fTransMat, myHMM.mTransMat) ;

	switch (myDistrType)
	{	case eNormalDistr :
		{	cUnivariateNormal *myLoi = (cUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{	cMultivariateNormal *myLoi = (cMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNbClasses, myLoi->mCov) ;
		}
		break ;
		case  eMixtUniNormalDistr :
		{	cMixtUnivariateNormal *myParam = (cMixtUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNbClasses, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNbClasses, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eMixtMultiNormalDistr :
		{	cMixtMultivariateNormal *myParam = (cMixtMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNbClasses, myNbMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNbClasses, myNbMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{	cDiscrete *myParam = (cDiscrete *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myParam->mProba) ;
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

cOTMatrix* myProbaCond = new cOTMatrix[myNbSample] ;

	for (register uint n = 0 ; n < myNbSample ; n++)
		myProbaCond[n].ReAlloc(myNbClasses, myT[n]) ;

		myHMM.mDistrParam->ComputeCondProba(myY, myNbSample, myProbaCond) ;

cBaumWelch myBaumWelch=cBaumWelch(myNbSample, myT, myNbClasses) ;
	myBaumWelch.ForwardBackward(myProbaCond, myHMM) ;
// On enlève le scale	
	for (register uint n = 0 ; n < myNbSample ; n++)
	{	for (register uint t = 0 ; t < myT[n] ; t++)
			for (register uint j = 0 ; j < myNbClasses ; j++)
				myBaumWelch.mAlpha[n][j][t] *= myBaumWelch.mRho[n][t] ;	
		double myScale = 1.0L ;
		for (register int t = (int)myT[n]-1 ; t >= 0 ; t--)
		{	myScale *= myBaumWelch.mRho[n][t] ;
			for (register uint j = 0 ; j < myNbClasses ; j++)
				myBaumWelch.mBeta[n][j][t] *= myScale ;	
		}
	}

	for (register uint n = 0 ; n < myNbSample ; n++)
	{	myProbaCond[n].Delete() ;
		myY[n].Delete() ;
	}

	delete [] myY ;

	delete [] myProbaCond ;

SEXP	myAux[6] ;
uint*	myLigne = new uint[myNbSample] ;


	for (register uint n = 0 ; n < myNbSample ; n++)
		myLigne[n] = myNbClasses ;

	myRUtil.SetListMatSexp(myBaumWelch.mAlpha, myNbSample,myAux[0]) ;
	myRUtil.SetListMatSexp(myBaumWelch.mBeta, myNbSample, myAux[1]) ;
	myRUtil.SetListMatSexp(myBaumWelch.mGamma, myNbSample, myAux[2]) ;
	myRUtil.SetListListMatSexp(myBaumWelch.mXsi, myNbSample, myT, myAux[3]) ;
	myRUtil.SetListVectSexp(myBaumWelch.mRho, myNbSample, myAux[4]) ;
	myRUtil.SetListValSexp(myBaumWelch.mLogVrais, myAux[5]) ;

	delete [] myLigne ;
	delete [] myT ;
SEXP myRes ;
	PROTECT(myRes = allocVector(VECSXP, 6)) ;
	for (register int i = 0 ; i < 6 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	myRUtil.EndProtect() ;

UNPROTECT(1) ;
	return(myRes) ;
}
END_EXTERN_C

BEG_EXTERN_C
DECL_DLL_EXPORT SEXP RLogforwardbackward	(	SEXP theHMM, 
												SEXP theYt
											)
{
distrDefinitionEnum		myDistrType	;
uint					myDimObs=1,
						myNbClasses,
						myNbProba=0,
						myNbMixt=0	;
cRUtil					myRUtil			;


SEXP myDistSEXP ;

	myRUtil.GetValSexp(theHMM, fDistr, myDistSEXP) ; // Loi de proba	
char myString[255] ;
char *myStr = (char *)myString ;
	myRUtil.GetValSexp(myDistSEXP, gType, myStr) ;
	myRUtil.GetValSexp(myDistSEXP, gNClasses, myNbClasses) ;
	if (strcmp(myStr, "NORMAL") == 0)
	{	myRUtil.GetValSexp(myDistSEXP, 2, myDimObs) ;
		if (myDimObs == 1)
			myDistrType = eNormalDistr ;
		else
			myDistrType = eMultiNormalDistr ;
	}
	else
	{	if (strcmp(myStr, "DISCRETE") == 0)
		{	myDistrType = eDiscreteDistr ;
			myRUtil.GetValSexp(myDistSEXP, 2, myNbProba) ;
		}
		else
		{	if (strcmp(myStr, "MIXTURE") == 0)
			{	myRUtil.GetValSexp(myDistSEXP, 2, myNbMixt) ;
				myRUtil.GetValSexp(myDistSEXP, 3, myDimObs) ;
				if (myDimObs == 1)
					myDistrType = eMixtUniNormalDistr ;
				else
					myDistrType = eMixtMultiNormalDistr ;
			}
		}
	}
uint	myNbSample = length(theYt) ;	
uint*	myT = new uint[myNbSample] ;

cOTVector* myY = new cOTVector[myNbSample] ;


	for (register uint n = 0 ; n < myNbSample ; n++)
	{	SEXP myAux ;
		myRUtil.GetValSexp(theYt, n, myAux) ;
		myT[n] = length(myAux) / myDimObs ;
		myY[n].ReAlloc(myT[n]*myDimObs) ;
		myY[n]= REAL(myAux) ;
	}

cHmm myHMM = cHmm(myDistrType, myNbClasses, myDimObs, myNbMixt, myNbProba) ;
	myRUtil.GetVectSexp(theHMM, fInitProba, myHMM.mInitProba) ;
	myRUtil.GetMatSexp(theHMM, fTransMat, myHMM.mTransMat) ;

	switch (myDistrType)
	{	case eNormalDistr :
		{	cUnivariateNormal *myLoi = (cUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetVectSexp(myDistSEXP, 3, myLoi->mMean) ;
			myRUtil.GetVectSexp(myDistSEXP, 4, myLoi->mVar) ;
		}
		break ;
		case eMultiNormalDistr :
		{	cMultivariateNormal *myLoi = (cMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myLoi->mMean) ; 
			myRUtil.GetListMatSexp(myDistSEXP, 4, myNbClasses, myLoi->mCov) ;
		}
		break ;
		case  eMixtUniNormalDistr :
		{	cMixtUnivariateNormal *myParam = (cMixtUnivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 4, myNbClasses, myParam->mMean) ;
			myRUtil.GetListVectSexp(myDistSEXP, 5, myNbClasses, myParam->mVar) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eMixtMultiNormalDistr :
		{	cMixtMultivariateNormal *myParam = (cMixtMultivariateNormal *)(myHMM.mDistrParam) ;
			myRUtil.GetListListVectSexp(myDistSEXP, 4, myNbClasses, myNbMixt, myParam->mMean) ;
			myRUtil.GetListListMatSexp(myDistSEXP, 5, myNbClasses, myNbMixt, myParam->mCov) ;
			myRUtil.GetListVectSexp(myDistSEXP, 6, myNbClasses, myParam->mp) ;
		}
		break ;

		case eDiscreteDistr :
		{	cDiscrete *myParam = (cDiscrete *)(myHMM.mDistrParam) ;
			myRUtil.GetListVectSexp(myDistSEXP, 3, myNbClasses, myParam->mProba) ;
		}
		break ;
		case eUnknownDistr :
		default :
		break ;
	}

cOTMatrix* myProbaCond = new cOTMatrix[myNbSample] ;

	for (register uint n = 0 ; n < myNbSample ; n++)
		myProbaCond[n].ReAlloc(myNbClasses, myT[n]) ;

		myHMM.mDistrParam->ComputeCondProba(myY, myNbSample, myProbaCond) ;

cLogBaumWelch myLogBaumWelch=cLogBaumWelch(myNbSample, myT, myNbClasses) ;
	myLogBaumWelch.LogForwardBackward(myProbaCond, myHMM) ;

	for (register uint n = 0 ; n < myNbSample ; n++)
	{	myProbaCond[n].Delete() ;
		myY[n].Delete() ;
	}

	delete [] myY ;

	delete [] myProbaCond ;

SEXP	myAux[6] ;
uint*	myLigne = new uint[myNbSample] ;


	for (register uint n = 0 ; n < myNbSample ; n++)
		myLigne[n] = myNbClasses ;

	myRUtil.SetListMatSexp(myLogBaumWelch.mLogAlpha, myNbSample,myAux[0]) ;
	myRUtil.SetListMatSexp(myLogBaumWelch.mLogBeta, myNbSample, myAux[1]) ;
	myRUtil.SetListMatSexp(myLogBaumWelch.mLogGamma, myNbSample, myAux[2]) ;
	myRUtil.SetListListMatSexp(myLogBaumWelch.mLogXsi, myNbSample, myT, myAux[3]) ;
	myRUtil.SetListVectSexp(myLogBaumWelch.mLogRho, myNbSample, myAux[4]) ;
	myRUtil.SetListValSexp(myLogBaumWelch.mLogVrais, myAux[5]) ;

	delete [] myLigne ;
	delete [] myT ;
SEXP myRes ;
	PROTECT(myRes = allocVector(VECSXP, 6)) ;
	for (register int i = 0 ; i < 6 ; i++)
		SET_VECTOR_ELT(myRes, i, myAux[i]) ;
	myRUtil.EndProtect() ;

UNPROTECT(1) ;
	return(myRes) ;
}
END_EXTERN_C
#endif //_RDLL_
