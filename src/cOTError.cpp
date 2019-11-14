/**************************************************************
 *** RHmm version 1.3.4                                      
 ***                                                         
 *** File: cOTError.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 ***                                                         
 *** Date: 2010/11/14                                      
 ***                                                         
 **************************************************************/

#include "cOTError.h"

cOTError::cOTError(char *theMess)
{
	if (theMess != (char *)NULL) 
		std::cout << theMess << std::endl ;
	exit(0) ;
}

