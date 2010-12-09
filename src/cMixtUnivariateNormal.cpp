/**************************************************************
 *** RHmm version 1.4.4                                     
 ***                                                         
 *** File: cMixtUnivariateNormal.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** Date: 2010/12/09                                     
 ***                                                         
 **************************************************************/

#include "cMixtUnivariateNormal.h"
static void MixtUnivariateNormalDensity(cOTVector& theY, uint theNMixt, cOTVector& theMean, cOTVector& theVar, cOTVector& thep, double* theDens)
{
cOTVector mySigma(theNMixt) ;
        for (register uint i = 0 ; i < theNMixt ; i++)
                mySigma[i] = sqrt(theVar[i]) ;

        for (register uint t = 0 ; t < theY.mSize ; t++)
        {       theDens[t] = 0.0 ;
                for (register uint i = 0 ; i < theNMixt ; i++)
                {
                double myCR = (theY[t] - theMean[i])/mySigma[i] ;
                        theDens[t] += thep[i]/(SQRT_TWO_PI*mySigma[i])*exp(-0.5*myCR*myCR) ;
                }
                if (theDens[t] < 1e-30)
                        theDens[t] = 1e-30 ;
        }
}

cMixtUnivariateNormal::cMixtUnivariateNormal(uint theNClass, uint theNMixt)
{       MESS_CREAT("cMixtUnivariateNormal")
        mvNClass = theNClass ;
        mvNMixt = theNMixt ;
        if ( (theNClass > 0) && (theNMixt > 0) )
        {       mMean = new cOTVector[theNClass] ;
                mVar = new cOTVector[theNClass] ;
                mp = new cOTVector[theNClass] ;
                for (register uint i = 0 ; i < mvNClass ; i++)
                {       mMean[i].ReAlloc(theNMixt) ;
                        mVar[i].ReAlloc(theNMixt) ;
                        mp[i].ReAlloc(theNMixt) ;
                }
        }
        else
        {       mMean = mVar = mp = NULL ;
                mvNClass = mvNMixt = 0 ;
        }
}

cMixtUnivariateNormal::~cMixtUnivariateNormal()
{       MESS_DESTR("cMixtUnivariateNormal")
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       mMean[i].Delete() ;
                mVar[i].Delete() ;
                mp[i].Delete() ;
        }
        mMean = mVar = mp = NULL ;
        mvNClass = mvNMixt = 0 ;
}

void cMixtUnivariateNormal::ComputeCondProba(cOTVector* theY, uint theNSample, cOTMatrix* theCondProba)
{
        for (register uint n = 0 ; n < theNSample ; n++)
                for (register uint i = 0 ; i < mvNClass ; i++)
                        MixtUnivariateNormalDensity(theY[n], mvNMixt, mMean[i], mVar[i], mp[i], theCondProba[n].mMat[i]) ;
}

void cMixtUnivariateNormal::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cOTMatrix* theCondProba)
{       
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       
        double mySumGammai = 0.0 ;
        register uint   n,
                                                t       ;
                for (n = 0 ; n < theInParam.mNSample ; n++)
                        for (t = 0 ; t < theInParam.mY[n].mSize  ; t++)
                                mySumGammai += theBaumWelch.mGamma[n][i][t] ;
                        
                for (register uint l = 0 ; l < mvNMixt ; l++)
                {                               
                        double myGammail ;
                        double mySumGammail = 0.0 ;
                        double myMoy = 0.0 ;
                        double myVar = 0.0 ;
                        for (n = 0 ; n < theInParam.mNSample ; n++)
                                for (t = 0 ; t < theInParam.mY[n].mSize ; t++)
                                {       
                                double  myStd = sqrt(mVar[i][l]),
                                                myCR = (theInParam.mY[n][t] - mMean[i][l])/myStd ;
                                        myGammail = theBaumWelch.mGamma[n][i][t] * mp[i][l] * 1/(SQRT_TWO_PI*myStd) * exp(-0.5*myCR*myCR) / theCondProba[n][i][t] ;
                                        mySumGammail += myGammail ; 
                                        myMoy += myGammail * theInParam.mY[n][t] ;
                                        myVar += myGammail * theInParam.mY[n][t] * theInParam.mY[n][t] ;
                                }
                        mp[i][l] = mySumGammail / mySumGammai ;
                        mMean[i][l] = myMoy/mySumGammail ;
                        mVar[i][l] = myVar/mySumGammail ;
                        mVar[i][l] -= mMean[i][l] * mMean[i][l] ;
                }
        }
}               


void cMixtUnivariateNormal::InitParameters(cBaumWelchInParam &theInParam)
{
#ifdef _RDLL_
        GetRNGstate();
#endif //_RDLL_
double  myMoy = 0, 
                myVar = 0,
                mystdev         ;
register uint s = 0 ;
        for (register uint n = 0 ; n < theInParam.mNSample ; n++)
                for (register uint t = 0 ; t < theInParam.mY[n].mSize  ; t++)
                        {       myMoy = ((double)s*myMoy + theInParam.mY[n][t])/(double)(s+1) ;
                                myVar = ((double)s*myVar + theInParam.mY[n][t]*theInParam.mY[n][t])/(double)(++s) ;
                        }
        myVar -= myMoy*myMoy ;
        mystdev = sqrt(myVar) ;

        for (register uint i = 0 ; i < mvNClass ; i++)
                {       double mySomme = 0.0 ;
                        register uint l ;
                        for (l = 0 ; l < mvNMixt ; l++)
                        {       mMean[i][l] =  -2*mystdev + myMoy + 2*mystdev * unif_rand() ;
                                mVar[i][l] = 0.5*myVar + 3*myVar * unif_rand() ;         
                                mp[i][l] = unif_rand() ;        
                                mySomme += mp[i][l] ;
                        }
                        for (l = 0 ; l < mvNMixt ; l++)
                                mp[i][l] /= mySomme ;
                }
                        
#ifdef _RDLL_
        PutRNGstate() ;
#endif //_RDLL_
}

void cMixtUnivariateNormal::CopyDistr(cDistribution* theSrc)
{
cMixtUnivariateNormal *mySrc ;
        mySrc = static_cast<cMixtUnivariateNormal *>(theSrc) ;
        mvNClass = mySrc->mvNClass ;
        mvNMixt = mySrc->mvNMixt ;
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       mMean[i] = mySrc->mMean[i] ;
                mVar[i] = mySrc->mVar[i] ;
                mp[i] = mySrc->mp[i] ;
        }
}

void cMixtUnivariateNormal::Print()
{
        Rprintf("Parameters\n") ;
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       Rprintf("State %d\n", i) ;
                for (register uint j = 0 ; j < mvNMixt ; j++)
                        Rprintf("\tm[%d]=%lf - s[%d]=%lf - p[%d]=%lf\n", j, mMean[i][j], 
                                j, sqrt(mVar[i][j]), j, mp[i][j]) ;
                Rprintf("\n") ;
        }
}


void cMixtUnivariateNormal::GetParam(uint theDeb, cOTVector& theParam)
{
register uint k = theDeb ;
        for (register uint n = 0 ; n < mvNClass ; n++)
                for (register uint p = 0 ; p < mvNMixt ; p++)
                {       theParam[k++] = mMean[n][p] ;
                        theParam[k++] = mVar[n][p] ;
                        if (p < mvNMixt-1)
                                theParam[k++] = mp[n][p] ;
                }       
}
void cMixtUnivariateNormal::SetParam(uint theDeb, cOTVector& theParam)
{
register uint k = theDeb ;
        for (register uint n = 0 ; n < mvNClass ; n++)
        {       mp[n][mvNMixt-1] = 1.0L ;
                for (register uint p = 0 ; p < mvNMixt ; p++)
                {       mMean[n][p] = theParam[k++] ;
                        mVar[n][p] = theParam[k++] ;
                        if (p < mvNMixt-1)
                        {       mp[n][p] = theParam[k++] ;
                                mp[n][mvNMixt-1] -= mp[n][p] ;
                        }
                }
        }
}

