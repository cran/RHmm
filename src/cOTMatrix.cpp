/***********************************************************
 * RHmm version 0.9.3                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/11/07                                        *
 *                                                         *
 ***********************************************************/
#include "cOTMatrix.h"


cOTVector& cOTVector::operator =(cOTMatrix& theMatrix)
{
	if (theMatrix.mNCol == 1)
	{
		if (mSize == 0)
		{	mVect = new double[theMatrix.mNRow] ;
			mSize = theMatrix.mNRow ;
		}
		else
		{	delete [] mVect ;
			mSize = theMatrix.mNRow ;
			mVect = new double[mSize] ;
		}
		for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] = theMatrix.mMat[i][0] ;
	}
	return *this ;
}


cOTMatrix& transpose(cOTVector &theVect)
{
cOTMatrix* myTranspose ;
	myTranspose = new cOTMatrix(1, theVect.mSize) ;
	for (register uint i=0 ; i < theVect.mSize ; i++)
			myTranspose->mMat[0][i] = theVect[i] ;
	return *myTranspose ;
}

cOTMatrix::cOTMatrix(uint theNRow, uint theNCol, double theVal)
{
	if ( (theNRow == 0) && (theNCol == 0) )
	{	mNRow = mNCol = 0 ;
		mMat = (double **)NULL ;
	}
	else
	{	if (theNRow > 0) 
		{	if ( (mMat = new double*[theNRow]) == NULL )
				throw cOTError("memory allocation problem") ;
			for (register uint i = 0 ; i < theNRow ; i++)
			{	if (theNCol > 0)
				{	if ( (mMat[i] = new double[theNCol]) == NULL )
						throw cOTError("memory allocation problem") ;
					for (register uint j = 0 ; j < theNCol ; j++)
						mMat[i][j] = theVal ;
				}
				else
					mMat[i] = (double *)NULL ;
			}
			mNRow = theNRow ;
			mNCol = theNCol ;
		}
		else
			throw cOTError("NRow must be strictly positive") ;
	}
}

cOTMatrix::~cOTMatrix()
{
	if (mNRow > 0)
	{	for (register uint i = 0 ; i < mNRow ; i++)
			delete [] mMat[i] ;
		delete [] mMat ;
		mMat = (double **)NULL ;
	}
	mNRow = mNCol = 0 ;
}

void cOTMatrix::Delete(void)
{
	for (register uint i = 0 ; i < mNRow ; i++)
		delete [] mMat[i] ;
	if (mMat != NULL)
		delete [] mMat ;
	mMat = (double **)NULL ;
	mNRow = mNCol = 0 ;
}

void cOTMatrix::ReAlloc(uint theNRow, uint theNCol, double theVal)
{
uint	i	;
	Delete() ;

	if ( (theNRow > 0) && (theNCol > 0) )
	{	if ( (mMat = new double*[theNRow]) == NULL)
			throw cOTError("memory allocation problem");
		mNRow = theNRow ;
		for (i = 0 ; i < mNRow ; i++)
			if ( (mMat[i] = new double[theNCol]) == NULL)
				throw cOTError("memory allocation problem") ;
			else
				for (register uint j = 0 ; j < theNCol ; j++)
					mMat[i][j] = theVal ;
		mNCol = theNCol ;
	}
}
double* & cOTMatrix::operator [](uint theNRow)
{
	if (theNRow < mNRow)
		return mMat[theNRow] ;
	else
		throw cOTError("bad index") ;
}


cOTMatrix& cOTMatrix::operator =(cOTMatrix& theSrcMatrix)
{
register uint	i,
				j	;
	
	if (mNRow == 0)
	{	if ( (mMat = new double*[theSrcMatrix.mNRow]) == NULL)
			throw cOTError("memory allocation problem") ;
		mNRow = theSrcMatrix.mNRow ;
		for (i = 0 ; i < mNRow ; i++)
			if ( (mMat[i] = new double[theSrcMatrix.mNCol]) == NULL)
				throw cOTError("memory allocation problem") ;
		mNCol = theSrcMatrix.mNCol ;
	}
	else
	{	if ( (mNRow != theSrcMatrix.mNRow) || (mNCol != theSrcMatrix.mNCol) )
		{	for (i = 0 ; i < mNRow ; i++)
				delete [] mMat[i]  ;
			delete [] mMat ;
			mNRow = theSrcMatrix.mNRow ;
			mNCol = theSrcMatrix.mNCol ;
			if ( (mMat = new double*[theSrcMatrix.mNRow]) == NULL)
				throw cOTError("memory allocation problem") ;
			for (i = 0 ; i < mNRow ; i++)
				if ( (mMat[i] = new double[theSrcMatrix.mNCol]) == NULL)
					throw cOTError("memory allocation problem") ;
		}
	}
	for (i = 0 ; i < mNRow ; i++)
		for (j = 0 ; j < mNCol ; j++)
			mMat[i][j] = theSrcMatrix.mMat[i][j] ;
	return *this ;
}
cOTMatrix& cOTMatrix::operator =(cOTVector& theVect)
{
register uint	i	;
	
	if (mNRow == 0)
	{	if ( (mMat = new double*[theVect.mSize]) == NULL)
			throw cOTError("memory allocation problem") ;
		mNRow = theVect.mSize ;
		for (i = 0 ; i < mNRow ; i++)
			if ( (mMat[i] = new double[1]) == NULL)
				throw cOTError("memory allocation problem") ;
		mNCol = 1 ;
	}
	else
	{	if ( (mNRow != theVect.mSize) || (mNCol != 1) )
		{	for (i = 0 ; i < mNRow ; i++)
				delete [] mMat[i]  ;
			delete [] mMat ;
			mNRow = theVect.mSize;
			mNCol = 1 ;
			if ( (mMat = new double*[mNRow]) == NULL)
				throw cOTError("memory allocation problem") ;
			for (i = 0 ; i < mNRow ; i++)
				if ( (mMat[i] = new double[1]) == NULL)
					throw cOTError("memory allocation problem") ;
		}
	}
	for (i = 0 ; i < mNRow ; i++)
		mMat[i][0] = theVect[i] ;
	return *this ;
}
cOTMatrix& cOTMatrix::operator =(double theVal)
{
	if ( (mNRow > 0) && (mNCol > 0) )
	{	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
				mMat[i][j] = theVal ;
	}
	return *this ;
}
cOTMatrix& cOTMatrix::operator +(cOTMatrix& theMatrix)
{	if ( (theMatrix.mNCol == mNCol) && (theMatrix.mNRow == mNRow) )
	{	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
				mMat[i][j] += theMatrix.mMat[i][j] ;
	}
	return *this ;
}
cOTMatrix& cOTMatrix::operator +=(cOTMatrix& theMatrix)
{	if ( (theMatrix.mNCol == mNCol) && (theMatrix.mNRow == mNRow) )
	{	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
				mMat[i][j] += theMatrix.mMat[i][j] ;
		return *this ;
	}
	else
		throw cOTError("wrong matrices size") ;
}

cOTMatrix& cOTMatrix::operator -(cOTMatrix& theMatrix)
{	if ( (theMatrix.mNCol == mNCol) && (theMatrix.mNRow == mNRow) )
	{	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
				mMat[i][j] -= theMatrix.mMat[i][j] ;
		return *this ;
	}
	else
		throw cOTError("wrong matrices size") ;
}
cOTMatrix& cOTMatrix::operator -=(cOTMatrix& theMatrix)
{	if ( (theMatrix.mNCol == mNCol) && (theMatrix.mNRow == mNRow) )
	{	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
				mMat[i][j] -= theMatrix.mMat[i][j] ;
		return *this ;
	}
	else
		throw cOTError("wrong matrices size") ;
}


cOTMatrix& operator *(cOTMatrix& theLeft, cOTMatrix &theRight)
{	
cOTMatrix *myRes = new cOTMatrix(theLeft.mNRow, theRight.mNCol) ;
	if ( (theLeft.mNCol == theRight.mNRow) )
	{	for (register uint i = 0 ; i < theLeft.mNRow ; i++)
			for (register uint j = 0 ; j < theRight.mNCol ; j++)
				for (register uint k = 0 ; k < theLeft.mNCol ; k++)
					myRes->mMat[i][j] += theLeft[i][k] * theRight[k][j] ;
		return *myRes ;
	}
	else
		throw cOTError("wrong matrices size") ;
}

cOTMatrix& cOTMatrix::operator *=(cOTMatrix& theMatrix)
{	
cOTMatrix myRes(mNRow, mNRow, 0.0L) ;
	if ( (theMatrix.mNCol == mNCol) && (theMatrix.mNRow == mNRow) && (mNRow == mNCol) )
	{	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
				for (register uint k = 0 ; k < mNRow ; k++)
					myRes[i][j] += mMat[i][k]*theMatrix.mMat[k][j] ;
		*this = myRes ;			
		return *this ;
	}
	else
		throw cOTError("wrong matrices size") ;
}

cOTMatrix& operator-(cOTMatrix& theRight)
{
cOTMatrix* myRes = new cOTMatrix(theRight.mNRow, theRight.mNCol) ;
	for (register uint i = 0 ; i < theRight.mNRow ; i++)
		for (register uint j = 0 ; j < theRight.mNCol ; j++)
			myRes->mMat[i][j]= -theRight.mMat[i][j] ;
	return *myRes ;
}



std::ostream& operator <<(std::ostream& theStream, cOTMatrix& theMat)
{
register uint	i,
				j	;
	for (i = 0 ; i < theMat.mNRow ; i++)
	{	for (j = 0 ; j < theMat.mNCol-1 ; j++)
			theStream << theMat[i][j] << "\t" ;
		theStream << theMat[i][j] << std::endl ;
	}
	return theStream ;
}
cOTVector& operator *(cOTMatrix& theLeft, cOTVector& theVect)
{
cOTVector* myVect=(cOTVector *)NULL ;
	if (theLeft.mNCol == theVect.mSize)
	{	myVect = new cOTVector(theLeft.mNRow) ;
		for (register uint i = 0 ; i < theLeft.mNRow ; i++)
			for (register uint k= 0 ; k < theLeft.mNCol ; k++)
				myVect->mVect[i] += theLeft[i][k] * theVect[k] ;
		return *myVect ;
	}
	else
		throw cOTError("wrong matrix or vector size") ;
}

cOTMatrix& operator *(cOTVector& theVect, cOTMatrix& theRight)
{
cOTMatrix* myMat=(cOTMatrix *)NULL ;
	if (theRight.mNRow == 1)
	{	myMat = new cOTMatrix(theVect.mSize, theRight.mNCol) ;
		for (register uint i = 0 ; i < theVect.mSize ; i++)
			for (register uint j = 0 ; j < theRight.mNCol ; j++)
				myMat->mMat[i][j] += theVect[i]*theRight[0][j] ;
		return *myMat ;
	}
	else
		throw cOTError("wrong matrix or vector size") ;	
}
cOTMatrix& operator *(cOTMatrix& theMatrix, double theLambda)
{	
cOTMatrix* myRes = new cOTMatrix(theMatrix.mNRow, theMatrix.mNCol) ;
	for (register uint i = 0 ; i < theMatrix.mNRow ; i++)
		for (register uint j=0 ; j< theMatrix.mNCol ; j++)
			myRes->mMat[i][j] = theLambda*theMatrix.mMat[i][j] ;
	return *myRes ;
}

cOTMatrix& operator *(double theLambda, cOTMatrix& theMatrix)
{
	return(theMatrix*theLambda) ;
}
cOTMatrix& cOTMatrix::operator *=(double theLambda)
{	
	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
					mMat[i][j] *= theLambda ;	
	return *this ;
}

cOTMatrix& cOTMatrix::operator /(double theLambda)
{	
	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
					mMat[i][j] /= theLambda ;	
	return *this ;
}


cOTMatrix& cOTMatrix::operator /=(double theLambda)
{	
	for (register uint i = 0 ; i < mNRow ; i++)
			for (register uint j = 0 ; j < mNCol ; j++)
					mMat[i][j] /= theLambda ;	
	return *this ;
}


cOTMatrix& transpose(cOTMatrix &theMatrix)
{
cOTMatrix* myTranspose ;
	myTranspose = new cOTMatrix(theMatrix.mNCol, theMatrix.mNRow) ;
	for (register uint i=0 ; i < theMatrix.mNRow ; i++)
		for (register uint j = 0 ; j < theMatrix.mNCol ; j++)
			myTranspose->mMat[j][i] = theMatrix.mMat[i][j] ;
	return *myTranspose ;
}
cOTMatrix& zeros(uint theN, uint theP)
{
cOTMatrix *myMat = new cOTMatrix(theN, theP) ;
	return *myMat ;
}
cOTMatrix& identity(uint theN)
{
cOTMatrix *myMat = new cOTMatrix(theN, theN) ;
	for (register uint i=0 ; i < theN ; i++)
		myMat->mMat[i][i] = 1.0L ;
	return *myMat ;
}

/*void svd(cOTMatrix &theMatrix, cOTMatrix &theU, cOTVector &theS, cOTMatrix &theV)
{
	theU = theMatrix ;
	theV = identity(theMatrix.mNCol) ;
	theS = 0.0 ;
int	myErr = svd(theMatrix.mMat, theMatrix.mNRow, theMatrix.mNCol, theU.mMat, theS.mVect, theV.mMat) ; 
	if (myErr > 0)
		throw cOTError("svd: no convergence after 500 iterations") ;
}
*/

cOTMatrix& diag(cOTVector &theVect)
{
cOTMatrix* myMat = new cOTMatrix(theVect.mSize, theVect.mSize) ;
	for (register uint i = 0 ; i < theVect.mSize ; i++)
		myMat->mMat[i][i] = theVect.mVect[i] ;

	return *myMat ;
}
/*
cOTMatrix& inv(cOTMatrix &theMatrix)
{
cOTMatrix	myU,
			myV	;
cOTVector	myS	;
	svd(theMatrix, myU, myS, myV) ;

	for (register uint i = 0 ; i < myS.mSize ; i++)
		if (fabs(myS[i]) < MIN_DBLE)
			throw cOTError("Non inversible matrix") ;
		else
			myS[i] = 1.0L/myS[i] ;
cOTMatrix	myMatS = diag(myS) ;
	
	return (transpose(myV) * myMatS * transpose(myU)) ;
}

*/
cOTMatrix& inv(cOTMatrix &theMatrix)
{
cOTMatrix	myInv = cOTMatrix(theMatrix.mNRow, theMatrix.mNCol) ;
double myDet ;

	LapackInvAndDet(theMatrix, myInv, myDet) ;
	if (fabs(myDet) < MIN_DBLE)
			throw cOTError("Non inversible matrix") ;
	return myInv ;
}
void LapackInvAndDet(cOTMatrix& theMatrix, cOTMatrix& theInvMatrix, double& theDet)
{
double *myAP = new double[theMatrix.mNCol*(theMatrix.mNCol + 1)/2],
		*myW = new double[theMatrix.mNCol],
		*myZ = new double[theMatrix.mNCol*theMatrix.mNCol],
		*myWork = new double[theMatrix.mNCol * 3] ;
int myInfo,
	myN = (int)(theMatrix.mNCol),
	myldz = (int)(theMatrix.mNCol) ;

	for (register int i = 0 ; i < theMatrix.mNCol ; i++)
		for (register int j = i ; j < theMatrix.mNCol ; j++)
			myAP[i+(j+1)*j/2]  = theMatrix[i][j] ;

	F77_NAME(dspev)("V", "U", &myN, myAP, myW, myZ, &myldz, myWork, &myInfo) ;

	if (myInfo != 0)
		throw cOTError("Non inversible matrix") ;
	theDet = 1.0L ;
cOTVector myInvEigenValue = cOTVector(theMatrix.mNCol) ;

cOTMatrix myEigenVector = cOTMatrix(theMatrix.mNCol, theMatrix.mNCol) ;
	for (register uint i = 0 ; i < theMatrix.mNCol ; i++)
	{	theDet *= myW[i] ;
		myInvEigenValue[i] = 1.0 /myW[i] ;
		for (register int j = 0 ; j < theMatrix.mNCol ; j++)
			myEigenVector[i][j] = myZ[i + j*theMatrix.mNCol] ;
	}
	theInvMatrix =  myEigenVector * diag(myInvEigenValue) * transpose(myEigenVector);
	
	delete myAP ;
	delete myW ;
	delete myZ ;
	delete myWork ;

}

