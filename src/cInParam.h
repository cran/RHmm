/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
 *                                                         *
 ***********************************************************/
#ifndef _CINPARAM_H_
#define _CINPARAM_H_
#include "hmm.h"


class cInParam
{	public :
		distrDefinitionEnum	mDistrType		; // Type de loi de proba
		uint				mNClass			; // Nombre de classes
		uint				mDimObs			; // Dimension des observations
		uint				mNMixt			; // Nombre de lois mélangées
		uint				mNProba			; // Nombre de proba discrètes
		uint				mNSample		; // Nombre d'échantillons
		cOTVector*			mY				; // Tableau mNSample x mT[i] des observations
	public :
		cInParam(uint theNSample, uint theDimObs, cOTVector* theY, distrDefinitionEnum theDistrType=eNormalDistr, uint theNClass=2, uint theNMixt=0, uint theNProba=0) ;
		virtual ~cInParam() ;
		cInParam & operator =(const cInParam &theSrc) ;
		virtual void Print(void) ;
} ;


#endif //_CINPARAM_H_
