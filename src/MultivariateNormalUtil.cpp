/***********************************************************
 * RHmm version 1.0.3                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/06/26                                        *
 *                                                         *
 ***********************************************************/
#include "MultivariateNormalUtil.h"

void SymetricInverseAndDet	(	cOTMatrix&	theMat,
								double&		theDet,
								cOTMatrix&	theInvMat
							)
{
	LapackInvAndDet(theMat, theInvMat, theDet) ;
}



void MultivariateNormalDensity	(	cOTVector&	thex,
									cOTVector&	theMu,
									cOTMatrix&	theInvCov,
									double		theDet,
									double*		theDens
								) 

{
register uint	i,
				j,
				t				;
double			myAux,
				myRapport		;

uint myDimObs = theMu.mSize ;	
	myRapport = pow(SQRT_TWO_PI*sqrt(theDet), (int)myDimObs) ;

uint myT = thex.mSize / myDimObs ;

	for ( t = 0 ; t < myT ; t++)
	{	myAux = 0.0 ;
		for (i = 0 ; i < myDimObs ; i++)
			for (j = 0 ; j < myDimObs ; j++)
				myAux += (thex[t+i*myT]-theMu[i]) * theInvCov[i][j] * (thex[t+j*myT]-theMu[j]) ;
		theDens[t] = exp(-0.5*myAux)/myRapport ;
	}
}
