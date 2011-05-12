/**************************************************************
 *** RHmm version 1.5.0
 ***                                                         
 *** File: cMixtMultivariateNormal.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

static void MixtMultivariateNormalDensity(cDVector& theY, uint theNMixt, cDVector* theMean, cDMatrix* theInvCov, cDVector& theDet, cDVector& thep, double* theDens)
{
uint myT = theY.mSize / theMean[0].mSize ;
double* myDens = new double[myT] ;
        
        for (register uint t = 0 ; t < myT ; t++)
                theDens[t] = 0.0 ;
        for (register uint j = 0 ; j < theNMixt ; j++)
        {       MultivariateNormalDensity(theY, theMean[j], theInvCov[j], theDet[j], myDens) ;
                for (register uint t = 0 ; t < myT ; t++)
                        theDens[t] += thep[j] * myDens[t] ;                     
        }
        for (register uint t = 0 ; t < myT ; t++)
                theDens[t] = MAX(theDens[t], 1e-30) ;

        delete[] myDens ;
}
cMixtMultivariateNormal::cMixtMultivariateNormal(uint theNClass, uint theNMixt, uint theDimObs)
{       MESS_CREAT("cMixtMultivariateNormal")
        mvNClass = theNClass ;
        mvNMixt = theNMixt ;
        mvDimObs = theDimObs ;
        if ( (theNClass > 0) && (theNMixt > 0) && (theDimObs > 0) )
        {       mMean = new cDVector*[theNClass] ;
                mCov = new cDMatrix*[theNClass] ;
                mp = new cDVector[theNClass] ;
                for (register uint i = 0 ; i < mvNClass ; i++)
                {       mMean[i] = new cDVector[theNMixt] ;
                        mCov[i] = new cDMatrix[theNMixt] ;
                        mp[i].ReAlloc(theNMixt) ;
                        for (register uint j = 0 ; j < theNMixt ; j++)
                        {       mMean[i][j].ReAlloc(theDimObs) ;
                                mCov[i][j].ReAlloc(theDimObs, theDimObs) ;
                        }
                }
        }
        else
        {       mMean = NULL ;
                mp = NULL ;
                mCov = NULL ;
                mvNClass = mvNMixt = mvDimObs = 0 ;
        }
}

cMixtMultivariateNormal::~cMixtMultivariateNormal()
{       MESS_DESTR("cMixtMultivariateNormal")
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       for (register uint j = 0 ; j < mvNMixt ; j++)
                {       mMean[i][j].Delete() ;
                        mCov[i][j].Delete() ;
                }
                mp[i].Delete() ;
        }
        delete [] mMean ;
        delete [] mCov ;
        delete [] mp ;
        mMean = NULL ;
        mCov = NULL ;
        mp = NULL ;
        mvNClass = mvNMixt = mvDimObs = 0 ;
}

void cMixtMultivariateNormal::ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)
{
cDMatrix* myInvCov = new cDMatrix[mvNMixt] ;
cDVector myDet = cDVector(mvNMixt) ;

        for (register uint i = 0 ; i < mvNMixt ; i++)
                myInvCov[i].ReAlloc(mvDimObs, mvDimObs) ;


        for (register uint i = 0 ; i < mvNClass ; i++)
        {       for (register uint j = 0 ; j < mvNMixt ; j++)
                        {
                        SymetricInverseAndDet(mCov[i][j], myDet[j], myInvCov[j]) ;
                        }
                for (register uint n = 0 ; n < theNSample ; n++)
                        MixtMultivariateNormalDensity(theY[n], mvNMixt, mMean[i], myInvCov, myDet, mp[i], theCondProba[n][i]) ;
        }
        for (register uint i = 0 ; i < mvNMixt ; i++)
                myInvCov[i].Delete() ;
        delete [] myInvCov ;
}

void cMixtMultivariateNormal::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba)
{       
cDMatrix*      myInvCov = new cDMatrix[mvNMixt];
double*         myDet = new double[mvNMixt] ;

        for (register uint i = 0 ; i < mvNMixt ; i++)
                myInvCov[i].ReAlloc(mvDimObs, mvDimObs) ;

        for (register uint i = 0 ; i < mvNClass ; i++)
        {       
        double mySumGammai = 0.0 ;
        register uint   n,
                                        t       ;
                
                for (n = 0 ; n < theInParam.mNSample ; n++)
                {       uint myT = theInParam.mY[n].mSize / mvDimObs ;
                        for (t = 0 ; t < myT  ; t++)
                                mySumGammai += theBaumWelch.mGamma[n][i][t] ;
                }

                for (register uint j = 0 ; j < mvNMixt ; j++)
                        SymetricInverseAndDet(mCov[i][j], myDet[j], myInvCov[j]) ;
                
        cDVector myMoy = cDVector(mvDimObs) ;
        cDMatrix myCov = cDMatrix(mvDimObs, mvDimObs) ;
                for (register uint l = 0 ; l < mvNMixt  ; l++)
                {                               
                        myMoy = 0.0 ;
                        myCov = 0.0 ;
                        double myGammail ;
                        double mySumGammail = 0.0 ;
                        for (n = 0 ; n < theInParam.mNSample ; n++)
                        {       
                        uint myT = theInParam.mY[n].mSize / mvDimObs ;
                        double* myDens = new double[myT] ; 
                                MultivariateNormalDensity(theInParam.mY[n], mMean[i][l], myInvCov[l], myDet[l], myDens) ;
                                for (t = 0 ; t < myT ; t++)
                                {       myGammail = theBaumWelch.mGamma[n][i][t] * mp[i][l] * myDens[t] / theCondProba[n][i][t] ;
                                        mySumGammail += myGammail ; 
                                        for (register uint k = 0 ; k < mvDimObs ; k++)
                                        {       myMoy[k] += myGammail * theInParam.mY[n][t+k*myT] ;
                                                for (register uint j = k ; j < mvDimObs ; j++)
                                                        myCov[k][j] += myGammail * theInParam.mY[n][t+k*myT] * theInParam.mY[n][t+j*myT] ;
                                        }
                                }
                                delete [] myDens ;
                        }
                        mp[i][l] = mySumGammail / mySumGammai ;
                        mMean[i][l] = myMoy/mySumGammail ;
                        for (register int m = 0 ; m < (int)mvDimObs-1 ; m++)
                                for (register int l = m+1 ; l < (int)mvDimObs ; l++)
                                        myCov[l][m] = myCov[m][l] ;
                        mCov[i][l] = myCov/mySumGammail ;
                        mCov[i][l] -= mMean[i][l] * Transpose(mMean[i][l]) ;
                }
        }
}               


void cMixtMultivariateNormal::InitParameters(cBaumWelchInParam &theInParam)
{
#ifdef _RDLL_
        GetRNGstate();
#endif //_RDLL_

cDVector       myMoy(mvDimObs),
                        myVar(mvDimObs),
                        myStd(mvDimObs) ;

double mys = 0.0 ;
        for (register uint n = 0 ; n < theInParam.mNSample ; n++)
        {
        uint myT = theInParam.mY[n].mSize / mvDimObs ;
                for (register uint t = 0 ; t < myT  ; t++)
                {       for (register uint i = 0 ; i < mvDimObs ; i++)
                        {       myMoy[i] = (mys*myMoy[i] + theInParam.mY[n][t+i*myT])/(mys+1) ;
                                myVar[i] = (mys*myVar[i] + theInParam.mY[n][t+i*myT]*theInParam.mY[n][t+i*myT])/(mys+1) ;
                        }
                        mys++ ;
                }
        }

        for (register uint i = 0 ; i < mvDimObs ; i++)
        {       myVar[i] -= myMoy[i]*myMoy[i] ;
                myStd[i] = sqrt(myVar[i]) ;
        }

        for (register uint i = 0 ; i < mvNClass ; i++)
        {       
        double mySomme = 0.0 ;
        register uint l ;
                for (l = 0 ; l < mvNMixt ; l++)
                {
                                /* FIXME:
                                 * Zeros returns a global reference, but mCov[i][l] is a instance, so this
                                 * is fine, but weird of course.
                                 */
                                mCov[i][l] = Zeros(mCov[i][l].mNRow,mCov[i][l].mNCol);

                                for (register uint k = 0 ; k < mvDimObs ; k++)
                        {       mMean[i][l][k] =  -2*myStd[k] + myMoy[k] + 2*myStd[k] * unif_rand() ;
                                mCov[i][l][k][k] = 0.5*myVar[k] + 3*myVar[k] * unif_rand() ;     
                        }
                        mp[i][l] = unif_rand() ;        
                        mySomme += mp[i][l] ;
                }
                for (l = 0 ; l < mvNMixt ; l++)
                        mp[i][l] /= mySomme ;
        }

                        
#ifndef _RDLL_
        PutRNGstate() ;
#endif //_RDLL_
}

void cMixtMultivariateNormal::CopyDistr(cDistribution* theSrc)
{
cMixtMultivariateNormal *mySrc ;
        mySrc = static_cast<cMixtMultivariateNormal *>(theSrc) ;
        mvNClass = mySrc->mvNClass ;
        mvDimObs = mySrc->mvDimObs ;
        mvNMixt = mySrc->mvNMixt ;
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       for (register uint l = 0 ; l < mvNMixt ; l++)
                {       mMean[i][l] = mySrc->mMean[i][l] ;
                        mCov[i][l] = mySrc->mCov[i][l] ;
                }
                mp[i] = mySrc->mp[i] ;
        }
}

void cMixtMultivariateNormal::Print()
{
        Rprintf("Parameters\n") ;
        for (register uint i = 0 ; i < mvNClass ; i++)
        {       Rprintf("State %d\n", i) ;
                for (register uint j = 0 ; j < mvNMixt ; j++)
                {       Rprintf("p[%d]=%lf\nEsp[%d]\t\tMatCov[%d]\n", j, mp[i][j], j, j) ;
                        for (register uint k = 0 ; k < mvDimObs ; k++)
                        {       Rprintf("%lf\t", mMean[i][j][k]) ;
                                for (register uint l = 0 ; l < mvDimObs ; l++)
                                        Rprintf("\t%lf", mCov[i][j][k][l]) ;
                                Rprintf("\n") ;
                        }
                }
                Rprintf("\n") ;
        }
}


void cMixtMultivariateNormal::GetParam(uint theDeb, cDVector& theParam)
{
register uint k = theDeb ;
        for (register uint n = 0 ; n < mvNClass ; n++)
        {       for (register uint p = 0 ; p < mvNMixt ; p++)
                {       for (register uint m = 0 ; m < mvDimObs ; m++)
                                theParam[k++] = mMean[n][p][m] ;
                        for (register uint m = 0 ; m < mvDimObs ; m++)
                                for (register uint n = m ; n < mvDimObs ; n++)
                                        theParam[k++] = mCov[n][p][m][n] ;
                        if (p < mvNMixt-1)
                                theParam[k++] = mp[n][p] ;
                }
        }
}
void cMixtMultivariateNormal::SetParam(uint theDeb, cDVector& theParam)
{
register uint k = theDeb ;
        for (register uint n = 0 ; n < mvNClass ; n++)
        {       mp[n][mvNMixt-1] = 1.0 ;
                for (register uint p = 0 ; p < mvNMixt ; p++)
                {       for (register uint i = 0 ; i < mvDimObs ; i++)
                                mMean[n][p][i] = theParam[k++] ;
                        for (register uint i = 0 ; i < mvDimObs ; i++)
                                for (register uint j = i ; j < mvDimObs ; j++)
                                        mCov[n][p][i][j] = mCov[n][p][j][i] = theParam[k++] ;
                        if (p < mvNMixt-1)
                        {       mp[n][p] = theParam[k++] ;
                                mp[n][mvNMixt-1] -= mp[n][p] ;
                        }
                }
        }
}

