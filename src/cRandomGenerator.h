/**************************************************************
 *** RHmm version 1.4.2                                     
 ***                                                         
 *** File: cRandomGenerator.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/11/26                                     
 ***                                                         
 **************************************************************/

#ifndef _CRANDOMGENERATOR_H_
#define _CRANDOMGENERATOR_H_
#ifdef _GSL_
#include <gsl/gsl_rng.h>
        class cRandomGenerator
        {       public :
                        const gsl_rng_type      *mT     ;
                        gsl_rng                         *mr     ;
                public :
                        cRandomGenerator() ;
                        ~cRandomGenerator() ;
        } ;
#endif // _GSL_
#endif //_CRANDOMGENERATOR_H_
