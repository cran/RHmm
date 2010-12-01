/**************************************************************
 *** RHmm version 1.4.3                                     
 ***                                                         
 *** File: cOTError.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/01                                     
 ***                                                         
 **************************************************************/

#include "cOTError.h"

cOTError::cOTError(const char *theMess)
{
        if (theMess != (char *)NULL) 
                std::cout << theMess << std::endl ;
        exit(0) ;
}

