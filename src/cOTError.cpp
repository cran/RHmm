/***********************************************************
 * RHmm version 1.0.4                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2008/08/08                                        *
 *                                                         *
 ***********************************************************/
#include "coterror.h"

cOTError::cOTError(char *theMess)
{
	if (theMess != (char *)NULL) 
#ifndef _RDLL_
		cout << theMess << std::endl ;
#else
		Rprintf("%s\n", theMess) ;
#endif // _RDLL_
		exit(0) ;
}

