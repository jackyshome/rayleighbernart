/*
 * PffuncsRayleigh.cpp
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>

#include "smctc.hh"
#include "PffuncsRayleigh.h"

#include <boost/filesystem.hpp>

using namespace std;


    //dynamical model
    string flag_filter("Shannon");
    string saveDirSurrogate("/local/rdebroiz/tmp/");
    string saveDirTruePrior("/local/rdebroiz/tmp/");
    double domainLength(2.*M_PI);
    double alpha(1.);
    double timeStep(1.e-2);
    int nt(10);
    int nx(32);
    int ny(16);
    int nbSnapshots((int)nt/1);
    double RayNumber(1.e+1);//pow((double)(a_pos*a_pos+M_PI*M_PI),3.)/(a_pos*a_pos);


    //observation model variable
    int maxObs(200);
    int tStar(0);
    int tFreq(2);
    int SamplingRate(1);
    double NoiseVar(1.e-2);//(1.e-5);
    double sigmaNoise(0.e-2);
    CImg<double> TranfertFunction;

    // prior uniform in the following intervals:
    double xminTrue=-0.5,yminTrue=-0.5,zminTrue=0,aminTrue=M_PI/sqrt(2)-2*M_PI/2.,prminTrue=5.e+0;
    double xmaxTrue=0.5,ymaxTrue=0.5,zmaxTrue=1,amaxTrue=M_PI/sqrt(2)+2*M_PI/2.,prmaxTrue=5.e+0;

    // prior surrogate uniform in the following intervals:
    double xmin=-1,ymin=-1,zmin=0,amin=M_PI/sqrt(2)-2*M_PI/2,prmin=5.e+0;
    double xmax=1,ymax=1,zmax=2,amax=M_PI/sqrt(2)+2*M_PI/2,prmax=5.e+0;


///The function corresponding to the log likelihood at specified time and position (up to normalisation)
///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider
double logLikelihood(long lTime, cv_state & X)
{

    double res=0;
    CImg<double> currentObs,transfertFunction;
    CImgList<double> currentState(2);
    transfertFunction.assign((X.omObj).getTF());
    for (int k=0;k<maxObs;k++){
        currentObs.assign((X.omObj).getObs(lTime,k));
        if (currentObs.size()!=0){
            currentState.assign((X.dmObj).getCurrentfState());
            cimglist_for(currentState,l){currentState[l]=currentState[l].get_mul(transfertFunction.get_slice(0))*4*(X.dmObj).getNx()*(X.dmObj).getNy();}
            res+=(currentObs.get_slice(0)-(currentState.get_FFT(true))[0]).get_pow(2.).mean();
            currentState.assign((X.dmObj).getCurrentfState_th());
            cimglist_for(currentState,l){currentState[l]=currentState[l].get_mul(transfertFunction.get_slice(1))*4*(X.dmObj).getNx()*(X.dmObj).getNy();}
            res+=(currentObs.get_slice(1)-(currentState.get_FFT(true))[0]).get_pow(2.).mean();
        }
    }
    res*=-1.e0;
    cout<<"\n unormalized - log likelihood at time "<<lTime<<" = "<<-1*res<<", ";
    return res;
}

///A function to initialise particles
/// \param pRng A pointer to the random number generator which is to be used
smc::particle<cv_state> fInitialise(smc::rng *pRng)
{

    cv_state *xi= new cv_state[1];

    //dynamical model initialization
    xi->x_pos = pRng->Uniform(xmin,xmax);
    xi->y_pos = pRng->Uniform(ymin,ymax);
    xi->z_pos = pRng->Uniform(zmin,zmax);
    xi->a_pos = pRng->Uniform(amin,amax);
    xi->pr_pos = pRng->Uniform(prmin,prmax);


    xi->dmObj = DynamicalModel(nx,ny,nt,nbSnapshots,timeStep,xi->pr_pos,-1, domainLength,alpha,flag_filter,"");
    xi->dmObj.setSpatialInitialCondition(xi->dmObj.setLorenzInitialCondition(xi->x_pos,xi->y_pos,xi->z_pos,xi->a_pos,1.,xi->pr_pos));
    xi->dmObj.setRayleighNumber(RayNumber);
    cout<<"Parameters of the dynamical model:"<<endl;
    xi->dmObj.display(); cout<<endl;

    //observation model initialization
    xi->omObj= y->omObj;

    //path setting
    ostringstream oss1,oss2,oss3,oss4,oss5,oss6,oss7,oss8,oss9,oss10;
    oss1 <<nx;oss2 <<ny;oss3 <<nt;oss4<<timeStep;oss5<<xi->pr_pos;oss7<<xi->a_pos;oss8<<xi->x_pos;oss9<<xi->y_pos;oss10<<xi->z_pos;
    string mainDir= saveDirSurrogate + "_nx"+ oss1.str() + "_ny"+ oss2.str()+ "_nt"+ oss3.str()+ "_dt"+ oss4.str()+
                    "_Pr"+ oss5.str()+  "_a"+ oss7.str() + "_x"+ oss8.str()+ "_y"+ oss9.str()+ "_z"+ oss10.str()+"/";
    xi->dmObj.setSaveDir(mainDir);

    //creating new directories

    boost::filesystem::path saveDirSurrogatePath(saveDirSurrogate.c_str());
    boost::filesystem::create_directory(saveDirSurrogatePath);

    boost::filesystem::path mainDirPath(mainDir.c_str());
    boost::filesystem::create_directory(mainDirPath);

    double logp=logLikelihood(0,*xi);
    cout<<"1! "<<&(xi->dmObj)<<endl;
    return smc::particle<cv_state>(*xi,logp);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng)
{
    cv_state * cv_to = pFrom.GetValuePointer();
    //generate noise at time index "lTime"
    double amplitude=sigmaNoise*(int)(pRng->Uniform(0,2));// binary variable in {0,sigma}
    int *excitedFreq=new int[2];
    excitedFreq[0] = (int)(pRng->Uniform(1,(cv_to->dmObj).getNx()/20));// variable in {1,m_nx/10-1}
    excitedFreq[1] = 0;//(int)(pRng->Uniform(1,(cv_to->dmObj).getNy()/50));// variable in {1,m_ny/10-1}
    cout<<" ["<<amplitude<<","<<excitedFreq[0]<<","<<excitedFreq[1]<<"] ";
    //integrate one time step
    (cv_to->dmObj).setAdditionalNoise(amplitude,excitedFreq, lTime-1);
    (cv_to->dmObj).foward_oneStep(amplitude,excitedFreq,lTime);
    delete [] excitedFreq;
    //compute weihgts
    pFrom.AddToLogWeight(logLikelihood(lTime, *cv_to));
    cout<<" weight (unormalized) "<<pFrom.GetWeight()<<", ";
    //save state
    std::ostringstream oss;//,ossId;
    oss <<lTime;
    for (int tt=0;tt<lTime;tt++){
        //ossId<<(cv_to->dmObj).getAdditionalNoiseStd(tt);
        //ossId<<(cv_to->dmObj).getAdditionalNoiseFrq(tt)[0];
        //ossId<<(cv_to->dmObj).getAdditionalNoiseFrq(tt)[1];
        cout<<"\t ("<<(cv_to->dmObj).getAdditionalNoiseStd(tt)<<","<<(cv_to->dmObj).getAdditionalNoiseFrq(tt)[0]<<","<<(cv_to->dmObj).getAdditionalNoiseFrq(tt)[1]<<") ";
    }
    //(cv_to->dmObj).getCurrentState().save_cimg(const_cast<char*>(((cv_to->dmObj).getSaveDir()+"b" + oss.str() +"_"+ossId.str()+".cimg").c_str()));
    //(cv_to->dmObj).getCurrentState_th().save_cimg(const_cast<char*>(((cv_to->dmObj).getSaveDir()+"th" + oss.str() +"_"+ossId.str()+".cimg").c_str()));

}



void generate_observations(smc::rng *pRng){

    cv_state *xi= new cv_state[1];
    //dynamical model initialization
    DynamicalModel *dm=new DynamicalModel(nx,ny,nt,nbSnapshots,timeStep,-1.,-1, domainLength,alpha,flag_filter,"");
    xi->dmObj=*dm;
    //observation model initialization
    TranfertFunction.assign(2.*nx,2.*ny,2,1,1.).get_shared_slice(0).fill(0.);// only temperatures observations
    ObservationModel *om=new ObservationModel(&(xi->dmObj),maxObs, tStar, tFreq, SamplingRate, NoiseVar,TranfertFunction);
    xi->omObj=*om;
    cout<<"Parameters of the observation model:"<<endl;
    xi->omObj.display();cout<<endl;
    double x_pos, y_pos, z_pos, a_pos ,pr_pos ;

    for (int k=0;k<maxObs;k++){

        //dynamical model initialization
        x_pos = pRng->Uniform(xminTrue,xmaxTrue);
        y_pos = pRng->Uniform(yminTrue,ymaxTrue);
        z_pos = pRng->Uniform(zminTrue,zmaxTrue);
        a_pos = pRng->Uniform(aminTrue,amaxTrue);
        pr_pos = pRng->Uniform(prminTrue,prmaxTrue);
        xi->dmObj.setSpatialInitialCondition(xi->dmObj.setLorenzInitialCondition(x_pos,y_pos,z_pos,a_pos,1.,pr_pos));
        xi->dmObj.setRayleighNumber(RayNumber);
        cout<<"Parameters of the dynamical model:"<<endl;
        xi->dmObj.display(); cout<<endl;

        //path setting
        ostringstream oss1,oss2,oss3,oss4,oss5,oss6,oss7,oss8,oss9,oss10;
        oss1 <<nx;oss2 <<ny;oss3 <<nt;oss4<<timeStep;oss5<<pr_pos;oss7<<a_pos;oss8<<x_pos;oss9<<y_pos;oss10<<z_pos;

        string mainDir= saveDirTruePrior + "_nx"+ oss1.str() + "_ny"+ oss2.str()+ "_nt"+ oss3.str()+ "_dt"+ oss4.str()+
                    "_Pr"+ oss5.str()+  "_a"+ oss7.str() + "_x"+ oss8.str()+ "_y"+ oss9.str()+ "_z"+ oss10.str()+"/";

        string obsRootDir= saveDirTruePrior + "_nx"+ oss1.str() + "_ny"+ oss2.str()+ "_nt"+ oss3.str()+ "_dt"+ oss4.str()+"/";
        xi->dmObj.setSaveDir(mainDir);

        //creating new directories

        boost::filesystem::path saveDirTruePriorPath(saveDirTruePrior.c_str());
        if(boost::filesystem::create_directory(saveDirTruePriorPath))
            std::cout << "youpee";

        boost::filesystem::create_directory(boost::filesystem::path(mainDir.c_str()));
        string truthDir="stateTruth";
        boost::filesystem::create_directory(boost::filesystem::path((mainDir + truthDir)).c_str());

        string obsDir="observation";
        boost::filesystem::create_directory(boost::filesystem::path(obsRootDir.c_str()));
        boost::filesystem::create_directory(boost::filesystem::path((obsRootDir+obsDir)).c_str());

        //observation generation
        xi->dmObj.foward(truthDir);
        xi->omObj.generateObs(truthDir,obsRootDir+obsDir,k, pr_pos);

    }
    y=new cv_obs[1];
    y->omObj=xi->omObj;

}
