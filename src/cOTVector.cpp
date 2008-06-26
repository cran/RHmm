/***********************************************************
 * RHmm version 1.0.3                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/06/26                                        *
 *                                                         *
 ***********************************************************/
#include "cOTVector.h"

cOTVector::cOTVector()
{
	mSize = 0 ;
	mVect = (double *)NULL ;
}

cOTVector::cOTVector(uint theSize, double theVal)
{	mSize=theSize ; 
	if (theSize > 0)
	{	if ( (mVect = new double[theSize]) == NULL)
			throw cOTError("memory allocation problem") ; 
		for (register uint i = 0 ; i < theSize ; i++)
			mVect[i] = theVal ;
	}
	else
		mVect = (double *)NULL ;
}
cOTVector::cOTVector(uint theSize, double* theVect)
{	mSize=theSize ; 
	if (theSize > 0)
	{	if ( (mVect = new double[theSize]) == NULL)
			throw cOTError("memory allocation problem") ; 
		for (register uint i = 0 ; i < theSize ; i++)
			mVect[i] = theVect[i] ;
	}
	else
		mVect = (double *)NULL ;
}
cOTVector::~cOTVector()
{	if (mSize > 0)
	{	mSize = 0 ; 
		delete [] mVect ;
		mVect = (double *)NULL ;
	}
	else
		mVect = (double *)NULL ;
}
void cOTVector::Delete(void)
{
	if (mSize > 0)
		delete [] mVect ;
	mSize = 0 ;
	mVect = (double *)NULL ;
}

void cOTVector::ReAlloc(uint theSize)
{	Delete() ;
	if (theSize > 0)
	{	if ( (mVect = new double[theSize]) == NULL)
			throw cOTError("memory allocation problem") ; 
		mSize = theSize ;
	}
}

void cOTVector::ReAlloc(uint theSize, double theVal)
{	Delete() ;
	if (theSize > 0)
	{	if ( (mVect = new double[theSize]) == NULL)
			throw cOTError("memory allocation problem") ; 
		for (register uint n = 0 ; n < theSize ; n++)
			mVect[n] = theVal ;
		mSize = theSize ;
	}
}
void cOTVector::ReAlloc(uint theSize, double* theVect)
{	Delete() ;
	if (theSize > 0)
	{	if ( (mVect = new double[theSize]) == NULL)
			throw cOTError("memory allocation problem") ; 
		for (register uint n = 0 ; n < theSize ; n++)
			mVect[n] = theVect[n] ;
		mSize = theSize ;
	}
}

double& cOTVector::operator[](int theIndex)
{	if ( ((uint)theIndex >= 0) && ((uint)theIndex < mSize) )
		return(mVect[theIndex]) ;
	else
		throw cOTError("bad index") ;
}

cOTVector& cOTVector::operator =(cOTVector& theSrcVect)
{
	if (mSize == 0)
	{	if ( (mVect = new double[theSrcVect.mSize]) == NULL)
			throw cOTError("memory allocation problem") ;
		mSize = theSrcVect.mSize ;
	}
	else
	{	delete [] mVect ;
		mSize = theSrcVect.mSize ;
		mVect = new double[mSize] ;
	}
	for (register uint i = 0 ; i < mSize ; i++)
		mVect[i] = theSrcVect.mVect[i] ;
	return *this ;
}
cOTVector& cOTVector::operator =(double* theSrcVect)
{
	for (register uint i = 0 ; i < mSize ; i++)
		mVect[i] = theSrcVect[i] ;
	return *this ;
}

cOTVector& cOTVector::operator =(double theVal)
{
	if (mSize == 0)
	{	mSize = 1 ;
		mVect = new double[1] ;
		mVect[0] = theVal ;
	}
	else
	{	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] = theVal ;
	}
	return *this ;
}
cOTVector& cOTVector::operator +(cOTVector& theVect)
{
	if (mSize == theVect.mSize)
	{	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] += theVect.mVect[i] ;
		return *this ;
	}
	else
		throw cOTError("vectors must have the same size") ;

}
cOTVector& cOTVector::operator +(double theVal)
{
	if (mSize == 0)
	{	mSize = 1 ;
		mVect = new double[1] ;
		mVect[0] = 0.0L ;
	}

	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] += theVal ;
	
	return *this ;
}

cOTVector& cOTVector::operator +=(cOTVector& theSrcVect)
{
	if (mSize == theSrcVect.mSize)
	{	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] += theSrcVect.mVect[i] ;
		return *this ;
	}
	else
		throw cOTError("vectors must have the same size") ;
}

cOTVector& cOTVector::operator +=(double theVal)
{
	if (mSize == 0)
	{	mSize = 1 ;
		mVect = new double[1] ;
		mVect[0] = 0.0L ;
	}

	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] += theVal ;
	
	return *this ;
}


cOTVector& cOTVector::operator -(cOTVector& theVect)
{
	if (mSize == theVect.mSize)
	{	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] -= theVect.mVect[i] ;
		return *this ;
	}
	else
		throw cOTError("vectors must have the same size") ;
}
cOTVector& cOTVector::operator -(double theVal)
{
	if (mSize == 0)
	{	mSize = 1 ;
		mVect = new double[1] ;
		mVect[0] = 0.0L ;
	}

	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] -= theVal ;
	
	return *this ;
}
cOTVector& cOTVector::operator -=(cOTVector& theSrcVect)
{
	if (mSize == theSrcVect.mSize)
	{	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] -= theSrcVect.mVect[i] ;
		return *this ;
	}
	else
		throw cOTError("vectors must have the same size") ;

}
cOTVector& cOTVector::operator -=(double theVal)
{
	if (mSize == 0)
	{	mSize = 1 ;
		mVect = new double[1] ;
		mVect[0] = 0.0L ;
	}

	for (register uint i = 0 ; i < mSize ; i++)
			mVect[i] -= theVal ;
	
	return *this ;
}
cOTVector& cOTVector::operator *(double theLambda)
{	for (register uint i = 0 ; i < mSize ; i++)
		mVect[i] *= theLambda ;
	return *this ;
}
cOTVector& cOTVector::operator *=(double theLambda)
{	for (register uint i = 0 ; i < mSize ; i++)
		mVect[i] *= theLambda ;
	return *this ;
}
cOTVector& cOTVector::operator /(double theLambda)
{	for (register uint i = 0 ; i < mSize ; i++)
		mVect[i] /= theLambda ;
	return *this ;
}
cOTVector& cOTVector::operator /=(double theLambda)
{	for (register uint i = 0 ; i < mSize ; i++)
		mVect[i] /= theLambda ;
	return *this ;
}

bool operator ==(cOTVector& theVect1, cOTVector& theVect2)
{
	if (theVect1.mSize != theVect2.mSize)
		return(false) ;
	for (register uint i=0 ; i < theVect1.mSize ; i++)
		if (theVect1[i] != theVect2[i])
			return(false) ;
	return(true) ;
}

bool operator >(cOTVector& theVect1, cOTVector& theVect2)
{
	if (theVect1.mSize != theVect2.mSize)
		return(false) ;
	for (register uint i=0 ; i < theVect1.mSize ; i++)
		if (theVect1[i] <= theVect2[i])
			return(false) ;
	return(true) ;
}

bool operator >=(cOTVector& theVect1, cOTVector& theVect2)
{
	if (theVect1.mSize != theVect2.mSize)
		return(false) ;
	for (register uint i=0 ; i < theVect1.mSize ; i++)
		if (theVect1[i] < theVect2[i])
			return(false) ;
	return(true) ;
}
bool operator <(cOTVector& theVect1, cOTVector& theVect2)
{
	if (theVect1.mSize != theVect2.mSize)
		return(false) ;
	for (register uint i=0 ; i < theVect1.mSize ; i++)
		if (theVect1[i] >= theVect2[i])
			return(false) ;
	return(true) ;
}

bool operator <=(cOTVector& theVect1, cOTVector& theVect2)
{
	if (theVect1.mSize != theVect2.mSize)
		return(false) ;
	for (register uint i=0 ; i < theVect1.mSize ; i++)
		if (theVect1[i] > theVect2[i])
			return(false) ;
	return(true) ;
}



#ifndef _RDLL_
std::ostream& operator <<(std::ostream& theStream, cOTVector& theVect)
{
	for (register uint i = 0 ; i < theVect.mSize ; i++)
		theStream << theVect[i] << std::endl ;
	return theStream ;
}
#endif // _RDLL_

cOTVector& zeros(uint theN) 
{
cOTVector *myRes = new cOTVector(theN, 0.0L) ;
	return(*myRes) ;
}


cOTVector& copy_double(double *theVect, uint theSize)
{
cOTVector *myRes = new cOTVector(theSize, theVect) ;
	return *myRes ;
}