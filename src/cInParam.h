/**************************************************************
 *** RHmm version 1.4.3                                     
 ***                                                         
 *** File: cInParam.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/01                                     
 ***                                                         
 **************************************************************/

#ifndef _CINPARAM_H_
#define _CINPARAM_H_
#include "Hmm.h"


class cInParam
{       public :
                distrDefinitionEnum     mDistrType              ; // Type de loi de proba
                uint                            mNClass                 ; // Nombre de classes
                uint                            mDimObs                 ; // Dimension des observations
                uint                            mNMixt                  ; // Nombre de lois mélangées
                uint                            mNProba                 ; // Nombre de proba discrètes
                uint                            mNSample                ; // Nombre d'échantillons
                cOTVector*                      mY                              ; // Tableau mNSample x mT[i] des observations
        public :
                cInParam(uint theNSample, uint theDimObs, cOTVector* theY, distrDefinitionEnum theDistrType=eNormalDistr, uint theNClass=2, uint theNMixt=0, uint theNProba=0) ;
                virtual ~cInParam() ;
                cInParam & operator =(const cInParam &theSrc) ;
                virtual void Print(void) ;
} ;


#endif //_CINPARAM_H_
