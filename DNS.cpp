/*
 * DNS.cpp
 *
 *  Created on: 19 january 2016
 *      Author: pheas
 */

#include <iostream>
#include "DynamicalModel.h"
#include "ObservationModel.h"
#include "CImg.h"
using namespace std;
using namespace cimg_library;


int main () {


	string flag_filter("Shannon");
	string saveDir("./output/");
	double domainLength(2.*M_PI);
	double alpha(1.);
	int nx(128);
	int ny(64);
	int nt(1000/*10000*/);
	int nbSnapshots((int)nt/4);
	double timeStep(1.e-4);
	double diffCoef(1.e+1);
	double a(M_PI/sqrt(2));
	//double RaNumber(27*pow(M_PI,4.)/4.); // idem que RaNumber=pow(a*a+M_PI*M_PI,3.)/(a*a)
	double RaNumber(pow(a*a+M_PI*M_PI,3.)/(a*a)*40);


	//observation model variable
	int tStar(0);
	int tFreq(2);
	int SamplingRate(2);
	double NoiseVar(1.e-5);
	CImg<double> TranfertFunction;
	TranfertFunction.assign(2.*nx,2.*ny,2,1,1.).get_shared_plane(0).fill(0.);// only temperatures observations

	//path setting
	system(const_cast<char*>(("mkdir "+ saveDir).c_str()));
	ostringstream oss1,oss2,oss3,oss4,oss5,oss6,oss7,oss8;
	oss1 <<nx;oss2 <<ny;oss3 <<nt;oss4<<timeStep;oss5<<diffCoef;oss6<<domainLength;oss7<<RaNumber;oss8<<a;
	string mainDir= saveDir + "nx"+ oss1.str() + "_ny"+ oss2.str()+ "_nt"+ oss3.str()+ "_dt"+ oss4.str()+
			"_Pr"+ oss5.str()+ "_T"+ oss6.str()+ "_Ra"+ oss7.str()+ "_a"+ oss8.str() +"/";
 	system(const_cast<char*>(("mkdir "+ mainDir).c_str()));
	string truthDir="stateTruth";
	system(const_cast<char*>(("mkdir "+ mainDir+truthDir).c_str()));
	string obsDir="observation";
	system(const_cast<char*>(("mkdir "+ mainDir+obsDir).c_str()));
	string tmpDir="stateTmp";
	system(const_cast<char*>(("mkdir "+ mainDir+tmpDir).c_str()));
	string estimDir="stateEstim";
	system(const_cast<char*>(("mkdir "+ mainDir+estimDir).c_str()));
	string paramDir="parameters";
	system(const_cast<char*>(("mkdir "+ mainDir+paramDir).c_str()));


	DynamicalModel dm(nx,ny,nt,nbSnapshots,timeStep,diffCoef,RaNumber, domainLength,alpha,flag_filter,mainDir);
	ObservationModel om(&dm, 1,tStar, tFreq, SamplingRate, NoiseVar,TranfertFunction);
	cout<<"Parameters of the observation model:"<<endl;
	om.display();cout<<endl;

	//-launch time
	unsigned long clock;
	clock = cimg::time();

		 //-initialization
		 dm.setSpatialInitialCondition(dm.setLorenzInitialCondition(0.,1.,0.,a,40.,diffCoef));
		 cout<<"Parameters of the dynamical model:"<<endl;
		 dm.display(); cout<<endl;

		 dm.foward(truthDir);
		 //dm.fowardROMLorenz63(truthDir);

		 om.generateObs(truthDir,obsDir,0);
		/**/

	cout << "Total execution time : " << (cimg::time()-clock)/1000. << " s." << endl;
	cout << "****************************************************" << endl << endl;



    return 0;
}


