/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: MultivariateNormalUtil.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

void SymetricInverseAndDet(cDMatrix& theMat, double& theDet, cDMatrix& theInvMat)
{
        LapackInvAndDet(theMat, theInvMat, theDet) ;
}



void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, double*  theDens)
{
register uint i, j, t ;
double  myAux, myRapport ;

uint myDimObs = theMu.mSize ;   
        myRapport = pow(SQRT_TWO_PI, (int)myDimObs)*sqrt(theDet) ;

uint myT = thex.mSize / myDimObs ;

        for ( t = 0 ; t < myT ; t++)
        {       myAux = 0.0 ;
                for (i = 0 ; i < myDimObs ; i++)
                        for (j = 0 ; j < myDimObs ; j++)
                                myAux += (thex[t+i*myT]-theMu[i]) * theInvCov[i][j] * (thex[t+j*myT]-theMu[j]) ;
                theDens[t] = exp(-0.5*myAux)/myRapport ;
        }
}

void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, cDVector&  theDens)
{
register uint i, j, t ;
double  myAux, myRapport ;

uint myDimObs = theMu.mSize ;   
        myRapport = pow(SQRT_TWO_PI, (int)myDimObs)*sqrt(theDet) ;

uint myT = thex.mSize / myDimObs ;

        for ( t = 0 ; t < myT ; t++)
        {       myAux = 0.0 ;
                for (i = 0 ; i < myDimObs ; i++)
                        for (j = 0 ; j < myDimObs ; j++)
                                myAux += (thex[t+i*myT]-theMu[i]) * theInvCov[i][j] * (thex[t+j*myT]-theMu[j]) ;
                theDens[t] = exp(-0.5*myAux)/myRapport ;
        }
}
