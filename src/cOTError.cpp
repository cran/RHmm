/**************************************************************
 *** RHmm version 1.4.4                                     
 ***                                                         
 *** File: cOTError.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/09                                     
 ***                                                         
 **************************************************************/

#include "cOTError.h"

cOTError::cOTError(const char *theMess)
{
        if (theMess != (char *)NULL) 
                std::cout << theMess << std::endl ;
        std::exit(0) ;
}

