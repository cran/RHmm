/**************************************************************
 *** RHmm version 1.4.3                                     
 ***                                                         
 *** File: logprob.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/01                                     
 ***                                                         
 **************************************************************/

#include "logprob.h"

double eexp(const double theX)
{
	if (theX <= LOGZERO)
		return(0.0L) ;
	else
		return(exp(theX)) ;
}

double eln(const double theX)
{
	if (theX > 0.0L)
		return(log(theX)) ;
	else
		return(LOGZERO) ;
}

double elnsum1(const double theX, const double theY)
{	
double	myeLnX = eln(theX),
		myeLnY = eln(theY) ;
	
	if ( (myeLnX <= LOGZERO) || (myeLnY <= LOGZERO) )
	{	if (myeLnX <= LOGZERO)
			return(myeLnY) ;
		else
			return(myeLnX) ;
	}
	else
	{	if (myeLnX > myeLnY) 
			return(myeLnX + eln(1.0L+exp(myeLnY-myeLnX))) ;
		else
			return(myeLnY + eln(1.0L+exp(myeLnX-myeLnY))) ;
	}
}

double elnsum(const double theeLnX, const double theeLnY)
{	
// elnsum(eln(x), eln(y)) = eln(x+y) pour x, y > LOGZERO
// elnsum(LOGZERO, eln(y)) = eln(y)
// elnsum(eln(x), LOGZERO) = eln(x)
double	myeLnX = MAX(theeLnX, theeLnY),
		myeLnY = MIN(theeLnX, theeLnY) ;
	
	if (myeLnY <= LOGZERO)
		return(myeLnX) ;
	else
		return(myeLnX + eln(1.0L+exp(myeLnY-myeLnX))) ;
}


double elnproduct1(const double theX, const double theY)
{
double	myeLnX = eln(theX),
		myeLnY = eln(theY) ;

	if ( (myeLnX <= LOGZERO) || (myeLnY <= LOGZERO) )
		return(LOGZERO) ;
	else
		return(myeLnX + myeLnY) ;
}

double elnproduct(const double theeLnX, const double theeLnY)
// elnproduct(eln(x), eln(y)) = eln(x) + eln(y) pour x, y > 0
// elnproduct(LOGZERO, eln(y)) = elnproduct(eln(x), LOGZERO) = LOGZERO
{
	if ( (theeLnX <= LOGZERO) || (theeLnY <= LOGZERO) )
		return(LOGZERO) ;
	else
		return(theeLnX + theeLnY) ;
}
