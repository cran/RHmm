/***********************************************************
 * RHmm version 0.9.3                                      *
 *                                                         *
 *                                                         *
 * Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> *
 *                                                         *
 * Date: 2007/11/07                                        *
 *                                                         *
 ***********************************************************/
#ifndef _RHMM_H_
#define _RHMM_H_

typedef enum ParamHMMEnum
{	eNClasses=0,
	eObsDim,
	eNMixt,
	eNProba,
	eNoHmm,
	eDistrType
}ParamHMMEnum ;

typedef enum ParamAlgoBWEnum
{	eInitType=0,
	eNMaxIter,
	eTol,
	eVerbose,
	eNInitIter,
	eNMaxIterinit,
	eInitPoint
}ParamAlgoBWEnum ;

typedef enum HMMEnum
{	fInitProba=0,
	fTransMat,
	fDistr
}HMMEnum ;

typedef enum DistEnum
{	gType=0, 
	gNClasses
}DistEnum ;

#ifdef WIN32
	#define DECL_DLL_EXPORT __declspec(dllexport) 
#else
	#define DECL_DLL_EXPORT 
#endif // WIN32

#endif //_RHMM_H_

