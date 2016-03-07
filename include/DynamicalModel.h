/*
 * DynamicalModel.h
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#ifndef DYNAMICALMODEL_H_
#define DYNAMICALMODEL_H_



#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <time.h>
#include "CImg.h"
using namespace std;
using namespace cimg_library;

/*extern "C"{
#include "cblas.h"
#include "clapack.h"
}*/

class DynamicalModel {
public:
	DynamicalModel();
	DynamicalModel(int nx,
			int ny,
			int nt,
			int nbSnapshots,
			double timeStep,
			double diffCoef,
			double RaNumber,
			double domainLength,
			double alpha,
			string flag_filter,
			string mainDir
	);
	virtual ~DynamicalModel();
	DynamicalModel(const DynamicalModel& other);// copy constructor
	void copySubroutine(const DynamicalModel &original);
	DynamicalModel& operator = (const DynamicalModel &original); // assignement operator

	int getNx(); 					// returns the spatial resolution along the x-axis
	int getNy();					// returns the spatial resolution along the y-axis
	int getNt();					// returns the number of time steps
	int getNs();					// returns the number of snapshots
	double getTimeStep();			// returns the time step
	string getMainDir();			// returns main directory
	const CImgList<double>& getCurrentfState(); // returns fft b at current time
	const CImgList<double>& getInitialfState(); // returns fft th at initial time
	const CImg<double>& getCurrentState();// returns  b at current time
	const CImg<double>& getCurrentState_th();// returns th at current time
	const CImgList<double>& getCurrentfState_th(); // returns th at current time
	const CImgList<double>& getInitialfState_th(); // returns b at current time
	const CImgList<double>& getSnapshots(); // returns snapshots list
	int * getSnapshotIndices(); 				// returns snapshots temporal indices


	void setRayleighNumber(double Ra);
	void setPrandtlNumber(double Pr);
	const CImg<double>& setLorenzInitialCondition(double x0,double y0,double z0,double a63, double r63, double pr63); // set the value at time t0 of lorenz63 ROM
    void setSpatialInitialCondition(const CImg<double> initialCondition); // set the value of the state b and th at time t0
    void setFreqInitialCondition(const CImgList<double> initialCondition, const CImgList<double> initialCondition_th); // set the value of the state fftb and ffftth at time t0
	void setAdditionalNoise(double sigma, int *frequency, int timeIndex); // set the value of the additional noise at a given time index
	double getAdditionalNoiseStd(int timeIndex);
	int* getAdditionalNoiseFrq(int timeIndex);
	void setSaveDir(string saveDir);//  save directory
	string getSaveDir();
	void display();					// display the parameters of the dynamical model
	void foward(string saveDir);	// propagates the initial states m_b0 from t0 to m_nt*m_timeStep and saves the results in the directory "saveDir"
	void foward(string saveDir,const CImgList<double>&p); //POD-galerkin on subspace spanned by columns of p (no computational complexity reduction implemented here!)
	const CImgList<double>  & foward_oneStep(); 		// forward the state m_b (and m_fftb) from current_time to current_time+1 and update nextState
	const CImgList<double>  & foward_oneStep(double std, int *freq, int t); // forward the state m_b (and m_fftb) from current_time to current_time+1 and update nextState + add noise
	void fowardROMLorenz63(string saveDir);// propagates the initial states m_b0 from t0 to m_nt*m_timeStep with Lorenz63 ROM and saves the results in the directory "saveDir"
	void subspaceProj(const CImgList<double>&p, CImg<double> & b,CImg<double> & th,CImgList<double> & fftb,CImgList<double> & fftth); //proj on subspace spanned by columns of p

private:
	// Parameters
	 int m_nx;					// number of positive frequencies in the x-spatial dimension (# columns)
	 int m_ny;					// number of positive frequencies in the y-spatial dimension (# columns)
	 int m_nt;					// number of time instants (excluding initial time, ie 0,..,m_nt)
	 int m_nbSnapshots;					// nb of snapshots to backup
	 double m_timeStep;		// time step
	double m_diffCoef;		// coefficient of diffusion (or Prandtl number)
	double m_RaNumber;		// Rayleigh number
	 double m_domainLength;	// resolution domain is [0,m_domainLength]x[0,m_domainLength]
	 double m_alpha; 			// type of dynamical model ('alpha=0.5' for SQG, 'alpha=1.0' for NS2D)
	double *m_sigmaNoise;			//  standard deviation of additional noise (at each time step)
	int **m_freqNoise;			//  frequency of additional noise (at each time step)
	int *m_idx_snap;			// temporal indices of snapshots
	 string m_flag_filter;		// type of anti-aliasing filtering ('viscosity', viscosity_factor, 'hyper-vicosity', 'Shannon')
	string m_mainDir;

	// Working variables
	CImg<double> m_filter1;			// Filter designed to stabilize the numerical scheme (see [])
	CImg<double> m_filter2;			// Filter designed to stabilize the numerical scheme (see [])
	CImg<double> m_filter3;			// Filter designed to stabilize the numerical scheme (see [])

	CImg<double> m_b0;   			// b(x,y) on the resolution grid at t0
	CImgList<double> m_fftb0;		// b(k_x,k_y) fourier transform of b(x,y) at t0
	CImg<double> m_th0;   			// th(x,y) on the resolution grid at t0
	CImgList<double> m_fftth0;		// th(k_x,k_y) fourier transform of th(x,y) at t0
	CImgList<double> m_snap;		// saved snapshots

	CImg<double> m_b;    	 // b(x,y) on the resolution grid at tk
	CImgList<double> m_fftb; // b(k_x,k_y) fourier transform of b(x,y) at tk
	CImg<double> m_th;    	 // th(x,y) on the resolution grid at tk
	CImgList<double> m_fftth; //th(k_x,k_y) fourier transform of th(x,y) at tk
	CImg<double> m_vx;  	 // x-component of the velocity field on the resolution grid at tk
	CImg<double> m_vy;  	// y-component of the velocity field on the resolution grid at tk
	CImg<double> m_dbx; 	// db(x,y)/dx sampled on the resolution grid at tk
	CImg<double> m_dby; 	// db(x,y)/dy sampled on the resolution grid at tk
	CImg<double> m_dthx; 	// dth(x,y)/dx sampled on the resolution grid at tk
	CImg<double> m_dthy; 	// dth(x,y)/dy sampled on the resolution grid at tk

	CImgList<double> m_fftadv_term1,m_fftadv_termTh1;        // terme advection 1 RK4
	CImgList<double> m_fftadv_term2,m_fftadv_termTh2;        // terme advection 2 RK4
	CImgList<double> m_fftadv_term3,m_fftadv_termTh3;        // terme advection 3 RK4
	CImgList<double> m_fftadv_term4,m_fftadv_termTh4;        // terme advection 4 RK4

	CImg<double> m_vx2;
	CImg<double> m_vy2;
	CImg<double> m_dbx2,m_dthx2;
	CImg<double> m_dby2,m_dthy2;

	CImg<double> m_vx3;
	CImg<double> m_vy3;
	CImg<double> m_dbx3,m_dthx3;
	CImg<double> m_dby3,m_dthy3;

	CImg<double> m_vx4;
	CImg<double> m_vy4;
	CImg<double> m_dbx4,m_dthx4;
	CImg<double> m_dby4,m_dthy4;

	CImg<double> m_laplacien_freq;  // 2*m_my-2*m_mx matrix containing the Laplacian operator in the frequency domain
	CImg<double> m_laplacienInv_freq; // 2*m_my-2*m_mx matrix containing the inverse Laplacian operator in the frequency domain
	CImg<double> m_laplacienFrac_freq;  // 2*m_my-2*m_mx matrix containing the minus Laplacian fractional operator in the frequency domain
	CImg<double> m_fx,m_fy;			// 2*n_mx and 2*m_my vectors containing the frequencies resolved by the numerical scheme (ordered as in the fft transform)
	CImgList<double> m_ddx,m_ddy;		// 2*m_my-2*m_mx matrix containing the derivative operator (resp. in the x and y direction) in the frequency domain

	CImgList<double> m_fftvTmp;
	CImgList<double> m_fftvTmp2;

	CImg<double> m_tmpL63;

	CImg<double> m_ButterflyTraj; // trajectory of Lorenz63 ROM
	double *m_XYZ,*m_XYZ0; // current and intial state of Lorenz63 ROM
	double m_a63,m_r63,m_b63; // parameter of Lorenz63 ROM
	CImg<double> m_filter1_Lorenz63;

	/*Methods*/
	void m_allocate(); 						   // allocate memory to the private variables
	void m_set_filters(); 					   // set the values of m_fx,m_fy,m_laplacien_freq,m_ddx,m_ddy,m_filter1,m_filter2,m_filter3
	void m_foward_CrankNicolson_RungeKutta4(); // forward the state m_b (and m_fftb) from current_time to current_time+m_timeStep
	void m_fowardROMLorenz63_RungeKutta4(); // forward X,Y,Z  from current_time to current_time+m_timeStep
	const CImg<double>&  m_highDimStateFromROMLorenz63(double xx,double yy, double zz); // reconstruction High dimensional state from ROM state
	void m_advtermEval(const CImgList<double> &fftb_k ,CImgList<double> &fftadv_term,const CImg<double> &vx,const CImg<double> &vy,CImg<double> &dbx,CImg<double> &dby); // compute the advection term needed in m_foward_CrankNicolson_RungeKutta4() and update the variables dbx,dby
	void m_RungeKuttaEval(const CImg<double> &tmp,CImg<double> &term2); // rayleighBernhart term in Lorenz63 Model.
	void m_velocityEval(const CImgList<double> &fftb_k ,CImg<double> &vx,CImg<double> &vy); // evaluation of the velocity from b
	void m_rayleighBernhartBTerm(const CImgList<double> &fftth_k ,CImgList<double> &fftadv_term);// add to the b advection the RayleighBernart term needed in m_foward_CrankNicolson_RungeKutta4()
	void m_rayleighBernhartThTerm(const CImgList<double> &fftb_k ,CImgList<double> &fftadv_term);// compute the th advection the RayleighBernart term needed in m_foward_CrankNicolson_RungeKutta4()
	void m_zeroCross(CImg<double> &inOut); // zero forcing
	void m_zeroCross(CImgList<double> &inOut); // zero forcing
	const CImg<double>& m_ifft_fft(const CImg<double> &in);// forcing periodic real and zero cross forcing

};

#endif /* DYNAMICALMODEL_H_ */
