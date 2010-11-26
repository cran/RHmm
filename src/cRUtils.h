/**************************************************************
 *** RHmm version 1.4.2                                     
 ***                                                         
 *** File: cRUtils.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/11/26                                     
 ***                                                         
 **************************************************************/

#ifndef _CRUTILS_H_
#define _CRUTILS_H_

#include "OTMathUtil.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <vector>

#ifndef uint
        typedef unsigned int uint ;
#endif // uint

class cRUtil
{       private :
                int     mvNbProtect     ;
        public :
                cRUtil(){mvNbProtect = 0 ;} ;
                void EndProtect(void){if (mvNbProtect > 0) {UNPROTECT(mvNbProtect); mvNbProtect = 0 ; }} ;
                ~cRUtil(){mvNbProtect = 0 ;};
                /*
                 *      Récupérer une seule valeur à partir d'une liste SEXP à la place n° theNum
                 */
                void GetValSexp(SEXP theSEXP, uint theNum, uint &theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, int &theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, double &theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, char* theVal) ;
                void GetValSexp(SEXP theSEXP, uint theNum, SEXP &theVal) ;
                /*
                *       Récupérer une vecteur à partir d'une liste SEXP à la place n° theNum
                */
                void GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, double* theVal) ;
                void GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, int* theVal) ;
                void GetVectSexp(SEXP theSEXP, uint theNum, uint theDim, uint* theVal) ;
                void GetVectSexp(SEXP theSEXP, uint theNum, cOTVector& theVal) ;
                /*
                *       Récupérer une matrice à partir d'une liste SEXP à la place n° theNum
                */
                void GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, int** theMat) ;
                void GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, uint** theMat) ;
                void GetMatSexp(SEXP theSEXP, uint theNum, uint theLigne, uint theCol, double** theMat) ;
                void GetMatSexp(SEXP theSEXP, uint theNum, cOTMatrix& theMat) ;
                void GetMatListSexp(SEXP theSEXP, uint theNum, std::vector<cOTMatrix> &theList);
                /*
                *       Récupérer l'ensemble des nombres dans une liste de nombres
                */
                void GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, int* theVal) ;
                void GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, uint* theVal) ;
                void GetListValSexp(SEXP theSEXP, uint theNum, uint theNElt, double* theVal) ;
                /*
                * Récuperer l'ensemble des vecteurs dans une liste de vecteur
                */
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, int** theVal) ;
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, uint** theVal) ;
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theDim, double** theVal) ;
                void GetListVectSexp(SEXP theSEXP, uint theNum, uint theNElt, cOTVector* theVal) ;
                /*
                *       Récupérer l'ensemble des matrices d'une liste de matrices
                */
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, int*** theVal) ;
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, uint*** theVal) ;
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, uint theLigne, uint theCol, double*** theVal) ;
                void GetListMatSexp(SEXP theSEXP, uint theNum, uint theNElt, cOTMatrix* theVal) ;
                
                /*
                 * Récupérer l'ensemble des vecteurs dans une liste de liste de vecteurs
                 */
                void GetListListVectSexp(SEXP theSEXP, uint theNum, uint theNList1, uint theNList2, cOTVector** theVect) ;

                /*
                 * Récupérer l'ensemble des vecteurs dans une liste de liste de matrices
                 */
                void GetListListMatSexp(SEXP theSEXP, uint theNum, uint theNList1, uint theNList2, cOTMatrix** theVect) ;

                /*
                *       Remplit une seule valeur dans un SEXP à la place n° theNum 
                */
                void set_val_sexp(int theVal, SEXP &theSEXP) ;
                void set_val_sexp(uint theVal, SEXP &theSEXP) ;
                void set_val_sexp(double theVal, SEXP &theSEXP) ;
                /*
                *       Remplit un vecteur de taille theDim dans un SEXP 
                */
                void SetVectSexp(int *theVect, uint theDim, SEXP &theSEXP) ;
                void SetVectSexp(uint *theVect, uint theDim, SEXP &theSEXP) ;
                void SetVectSexp(double *theVect, uint theDim, SEXP &theSEXP) ;
                void SetVectSexp(cOTVector& theVect, SEXP &theSEXP) ;
          /*
                *       Remplit une matrice de taille theLigne x theCol dans un SEXP
                */
                void SetMatSexp(int **theMat, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetMatSexp(uint **theMat, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetMatSexp(double **theMat, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetMatSexp(cOTMatrix& theMat, SEXP &theSEXP) ;
                /*
                * Remplit une liste de theDim Nombres dans un SEXP
                */
                void SetListValSexp(int* theVal, uint theDim, SEXP &theSEXP) ;
                void SetListValSexp(uint* theVal, uint theDim, SEXP &theSEXP) ;
                void SetListValSexp(double* theVal, uint theDim, SEXP &theSEXP) ;
                void SetListValSexp(cOTVector& theVal,  SEXP &theSEXP) ;
                /*
                * Remplit une liste de theNElt vecteur de taille theDim dans un SEXP
                */
                void SetListVectSexp(int** theVal, uint theNElt, uint theDim, SEXP &theSEXP) ;
                void SetListVectSexp(uint** theVal, uint theNElt, uint theDim, SEXP &theSEXP) ;
                void SetListVectSexp(double** theVal, uint theNElt, uint theDim, SEXP &theSEXP) ;
                /*
                * Remplit une liste de theNElt vecteurs de tailles différentes theDim[i] dans un SEXP
                */
                void SetListVectSexp(int** theVal, uint theNElt, uint* theDim, SEXP &theSEXP) ;
                void SetListVectSexp(uint** theVal, uint theNElt, uint *theDim, SEXP &theSEXP) ;
                void SetListVectSexp(double** theVal, uint theNElt, uint *theDim, SEXP &theSEXP) ;
                void SetListVectSexp(cOTVector* theVal, uint theNElt, SEXP &theSEXP) ;
                /*
                * Remplit une liste de theNElt matrice de taille theLigne x theCol dans un SEXP
                */
                void SetListMatSexp(int*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetListMatSexp(uint*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP) ;
                void SetListMatSexp(double*** theVal, uint theNElt, uint theLigne, uint theCol, SEXP &theSEXP) ;
                /*
                * Remplit une liste de theNElt matrice de tailles diffréntes theLigne[i] x theCol[i] dans un SEXP
                */
                void SetListMatSexp(int*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP) ;
                void SetListMatSexp(uint*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP) ;
                void SetListMatSexp(double*** theVal, uint theNElt, uint *theLigne, uint *theCol, SEXP &theSEXP) ;
                void SetListMatSexp(cOTMatrix* theVal, uint theNElt, SEXP &theSEXP) ;

                /*
                * Remplit une liste de theNList1 elements de listes de theNList2 elements de vecteurs dans un SEXP
                */
                void SetListListVectSexp(cOTVector** theVect, uint theNList1, uint theNList2, SEXP &theSEXP) ;

                /*
                * Remplit une liste de theNList1 elements de listes de theNList2 elements de matrices dans un SEXP
                */
                void SetListListMatSexp(cOTMatrix** theMat, uint theNList1, uint theNList2, SEXP &theSEXP) ;
                void SetListListMatSexp(cOTMatrix** theMat, uint theNList1, uint* theNList2, SEXP &theSEXP) ;

        } ;
#endif // _CRUTILS_H_


