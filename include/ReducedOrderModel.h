/*
 * ReducedOrderModel.h
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#ifndef ROM_H_
#define ROM_H_

#define cimg_use_lapack

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include "CImg.h"
#include "DynamicalModel.h"
#include "ObservationModel.h"


using namespace std;
using namespace cimg_library;


#ifdef __cplusplus
extern "C" int dgeev_(char *, char *, int *, double *,  int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
#else
extern int dgeev_(char *, char *, int *, double *,  int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
#endif


class ROM {

public:
	ROM(DynamicalModel *samples, double * weight, int Np, int dim, int mmse);
	virtual ~ROM();
	void pod();
	CImgList<double> get_pod();
	void pop();
	CImgList<double> get_pop();
	void display() const;
	void checkPopBasis() ;
	CImg<double> runPopRom(const CImg<double> &initCond, int Nt, bool saveHDTraj);
	CImg<double> runPopRom(const CImg<double> &initCond, int Nt, int Nsubdim, bool saveHDTraj);
	CImg<double> runPopRom_HighComplexity(const CImg<double> &initCond, int Nt);
	CImg<double> runPodRom(const CImg<double> &initCond,DynamicalModel &dmObj,double pr);
	CImg<double> runPodRom(const CImg<double> &initCond,DynamicalModel &dmObj, double pr, int Nsubdim);

private:
	int m_Np,m_Nt, m_Ns,m_Nx,m_Ny,m_Dim;
	int m_mmse;							// default 0, classical snapshots methods =1, truth= -1
	DynamicalModel * m_samples; 			// pointer to set of dynamical model object containing snapshots
	double * m_weight;						// pointer to particle weights
	CImgList<double> m_podBasis; 			// k-dimensional POD basis
	CImgList<double> m_popBasis; 			// n-dimensional POP bases
	CImg<double> m_popEigenValue;
	CImg<double> m_popRomBasis; 			// k-dimensional POP reduced basis
	CImgList<double> *m_snaps;				//  snapshots
	int *m_idx_snap;						// temporal indices of snapshots
	CImgList<double> m_snapsApprox;			//  snapshot approx
	CImgList<double> m_snapsRomApprox;  	//  snapshot approx low dim representation

	void m_normalize_qBasis();
	void m_normalize_snapshots();
	void m_FindEigenValAndVec_NonSym(const CImg<double> &A, CImg<double> &valEigen,CImgList<double> &vecEigen_);
	CImgList<double> m_eigen_ztz( int timeIndexMin,int Ns);
	CImgList<double> m_eigen_zzt_from_ztz(const CImgList<double>& eigen_ztz,  int timeIndexMin,int Ns);
	CImgList<double> m_eigen_zzt(int timeIndexMin,int Ns);
	CImg<double> m_fowardOneStep_popRom(const CImg<double> &buffer);
	CImg<double> m_fowardOneStep_pop(const CImg<double> &buffer);
	CImg<double> m_pop_projLowDim(const CImg<double> &buffer);
	CImg<double> m_pop_recHighDim(const CImg<double> &buffer);

};

#endif /* ROM_H_ */
