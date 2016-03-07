/*
 * ObservationModel.h
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#ifndef OBSERVATIONMODEL_H_
#define OBSERVATIONMODEL_H_

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include "CImg.h"
#include "DynamicalModel.h"

using namespace std;
using namespace cimg_library;

class ObservationModel {
public:
	ObservationModel();
	ObservationModel(DynamicalModel * dm,
					  const int nbObs,
					 int tStar,
					 int tFreq,
					 int SamplingRate,
					 double NoiseVar,
					 CImg<double> TranfertFunction);
	virtual ~ObservationModel();
	ObservationModel(const ObservationModel& other);// copy constructor
	ObservationModel& operator = (const ObservationModel &original); // assignement operator
    void allocateCopy(const ObservationModel& original);
	void generateObs(string stateDir,string saveDir,int kthObs, double pr);
	void generateObs(string stateDir,string saveDir,int kthObs);
	const CImg<double> getObs(int l, int kthObs);
	const int getNbObs();
	double getPrandtlNumber(int kthObs);
	const CImg<double> getTruth(int l, int kthObs);
	const CImg<double>& getTF();
	void setTranfertFunction(const CImg<double> &TranfertFunction);
	void display();

private:
	 int m_nbObs; 				// number of observation sequence
	 int m_tStar; 				// starting observation time (between'1' and 'n_step+1')
	 int m_tFreq;				// Frequency of observation
	 int m_SamplingRate;			// subsampling rate
	 double m_NoiseVar;			// noise variance
	CImg<double> m_TranfertFunction;  // sensor transfert function
    //int ObsIdx[]
	CImgList<double> *m_obsList;  // list of observed images
	CImgList<double> *m_trueList;  // true images
	double *m_pr;					// prandtl parameters
	string flag; 					// sampling method ('grid','random')
	DynamicalModel * m_dm; 			// pointeur to dynamical model object
	int *m_idx_obs; 				//observation index vector
};

#endif /* OBSERVATIONMODEL_H_ */
