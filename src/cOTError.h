/**************************************************************
 *** RHmm version 1.4.2                                     
 ***                                                         
 *** File: cOTError.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/11/26                                     
 ***                                                         
 **************************************************************/

#ifndef _COTERROR_H_
#define _COTERROR_H_
#include <cstdlib>
#include <iostream>

#ifndef NULL
        #define NULL 0
#endif // NULL

#ifndef uint
        typedef unsigned int uint ;
#endif // uint

class cOTError
{
        public :
                cOTError(char *theMess) ;
} ;

#endif //_COTERROR_H_
