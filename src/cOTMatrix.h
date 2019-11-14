/**************************************************************
 *** RHmm version 1.3.4                                      
 ***                                                         
 *** File: cOTMatrix.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2010/11/14                                      
 ***                                                         
 **************************************************************/

#ifndef _COTMATRIX_H_
#define _COTMATRIX_H_
#ifndef MIN_DBLE
	#define MIN_DBLE 1e-16L
#endif //MIN_DBLE


#ifndef uint
	typedef unsigned int uint ;
#endif //uint

#include <iostream>
#include <math.h>
#include "R_ext/Lapack.h"
#include "cOTError.h"
#include "cOTVector.h"

class cOTMatrix
{
public :
	uint		mNRow	;
	uint		mNCol	;
	double**	mMat	;
public :
	cOTMatrix(uint theNRow=0, uint theNCol=0, double theVal = 0.0L) ;
	virtual ~cOTMatrix() ;
	void Delete(void) ;
	void ReAlloc(uint theNRow, uint theNCol, double theVal = 0.0L) ;
	double*& operator[](uint theNRow) ;
	cOTMatrix& operator =(cOTMatrix& theSrcMat) ;
	cOTMatrix& operator =(cOTVector& theVect) ;
	cOTMatrix& operator =(double theVal) ;
	cOTMatrix& operator +(cOTMatrix& theMatrix) ;
	cOTMatrix& operator +=(cOTMatrix& theMatrix) ;
	cOTMatrix& operator -(cOTMatrix& theMatrix) ;
	cOTMatrix& operator -=(cOTMatrix& theMatrix) ;
	friend cOTMatrix& operator -(cOTMatrix& theRight) ;
	friend cOTMatrix& operator *(cOTMatrix& theLeft, cOTMatrix &theRight) ;
	friend cOTVector& operator *(cOTMatrix& theLeft, cOTVector& theVect) ;
	friend cOTMatrix& operator *(cOTVector& theVect, cOTMatrix& theRight) ;
	friend cOTMatrix& operator *(cOTMatrix& theMat, double theLambda) ;
	friend cOTMatrix& operator *(double theLambda, cOTMatrix& theMat) ;
	cOTMatrix& operator *=(cOTMatrix& theRight) ;
	cOTMatrix& operator *=(double theLambda) ;
	cOTMatrix& operator /(double theLambda) ;
	cOTMatrix& operator /=(double theLambda) ;
	friend std::ostream& operator <<(std::ostream& theStream, cOTMatrix &theMat) ;
	friend cOTMatrix& transpose(cOTMatrix &theMatrix) ;
	friend cOTMatrix& zeros(uint theN, uint theP) ;
	friend cOTMatrix& identity(uint theN) ;
//	friend void svd(cOTMatrix &theMatrix, cOTMatrix &theU, cOTVector &theS, cOTMatrix &theV) ;
	friend cOTMatrix& diag(cOTVector &theVect) ;
	friend cOTMatrix& inv(cOTMatrix &theMatrix) ;
	friend void LapackInvAndDet(cOTMatrix &theMatrix, cOTMatrix &theInvMatrix, double& theDet) ;
} ;

#endif // _COTMATRIX_H_


	
