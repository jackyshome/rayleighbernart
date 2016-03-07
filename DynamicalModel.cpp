/*
 * DynamicalModel.cpp
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#include "DynamicalModel.h"

DynamicalModel::DynamicalModel(){
    cout<<" default const dyn! "<<this<<endl;


}
DynamicalModel::DynamicalModel(int nx,
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
							   ):m_nx(nx), m_ny(ny), m_nt(nt),  m_nbSnapshots(nbSnapshots), m_timeStep(timeStep),m_diffCoef(diffCoef), m_RaNumber(RaNumber), m_domainLength(domainLength),
							     m_alpha(alpha), m_flag_filter(flag_filter), m_mainDir(mainDir)
{
    cout<<" const dyn! "<<this<<endl;
    m_a63=-1.;
	m_r63=-1.;
	m_b63=-1.;
	// checks
	if (m_nx<=0){cout<<"Error: m_nx<=0"<<endl;}
	if (m_ny<=0){cout<<"Error: m_ny<=0"<<endl;}
	if (m_nt<=0){cout<<"Error: m_nt<=0"<<endl;}
	if (m_timeStep<=0){cout<<"Error: m_timeStep<=0"<<endl;}
	if (m_domainLength<=0){cout<<"Error: m_domainLength<=0"<<endl;}

    int tFreq=(int)((double)m_nt/(double)m_nbSnapshots+0.5);
    int count=0;
    for (int i=0;i<=m_nt;i+=tFreq){++count;}
    while (count>m_nbSnapshots){++m_nbSnapshots;}
    cout<<"nb snapshots="<<m_nbSnapshots;

	m_allocate();

    for (int i=0;i<=m_nt;i++){m_idx_snap[i]=0;}
    count=0;
    for (int i=0;i<=m_nt;i+=tFreq){m_idx_snap[i]=count;++count;}
	cout<<" snapshots sampling index: ";
	for (int i=0;i<=m_nt;i++){cout<<m_idx_snap[i]<<" ";}

	m_set_filters();

}

DynamicalModel::~DynamicalModel() {
//    cout<<" destructor of dyna model!! "<<this<<endl;
//    delete [] m_XYZ0;
//    delete [] m_XYZ;
//    delete [] m_sigmaNoise;
//    for (int i=0;i<m_nt;i++){ delete [] m_freqNoise[i];}
//    delete [] m_freqNoise;
//    delete [] m_idx_snap;
}

DynamicalModel::DynamicalModel(const DynamicalModel& original){
    cout<<"copy constructor of dyna model!! "<<this<<"="<<&original<<endl;
	m_nx=original.m_nx;
	m_ny=original.m_ny;
	m_nt=original.m_nt;
	m_nbSnapshots=original.m_nbSnapshots;
	m_allocate();
    copySubroutine(original);
}
DynamicalModel& DynamicalModel::operator = (const DynamicalModel &original){
    cout<<"assignement op. of dyna model!! "<<this<<"="<<&original<<endl;
	//cout<<"operator = of dyna model!!"<<endl;
	if (this != &original){// protect against invalid self-assignment

      //if (original.m_b0.is_empty()){cout<<"Warning: object is copied without being initialized."<<endl;}
      //else{
			m_nx=original.m_nx;
			m_ny=original.m_ny;
			m_nt=original.m_nt;
			m_nbSnapshots=original.m_nbSnapshots;
			m_allocate();
		    copySubroutine(original);
      //}

	}



    return *this;
}
void DynamicalModel::copySubroutine(const DynamicalModel &original){

            m_timeStep=original.m_timeStep;
            m_diffCoef=original.m_diffCoef;
            m_RaNumber=original.m_RaNumber;
            m_domainLength=original.m_domainLength;
            m_alpha=original.m_alpha;
            m_flag_filter=original.m_flag_filter;
            m_mainDir=original.m_mainDir;

            m_a63=original.m_a63;
            m_r63=original.m_r63;
            m_b63=original.m_b63;

            for (int i=0;i<3;i++){
                m_XYZ0[i]=original.m_XYZ0[i];
                m_XYZ[i]=original.m_XYZ[i];
            }
            for (int i=0;i<m_nt;i++){
                       m_sigmaNoise[i]=original.m_sigmaNoise[i];
                       m_freqNoise[i][0]=original.m_freqNoise[i][0];m_freqNoise[i][1]=original.m_freqNoise[i][1];
             }
            for (int i=0;i<=m_nt;i++){m_idx_snap[i]=original.m_idx_snap[i];}

            //non-shared variables
            m_snap=original.m_snap;
            m_b0=original.m_b0;m_th0=original.m_th0;
            m_fftb0=original.m_fftb0;m_fftth0=original.m_fftth0;
            m_b=original.m_b;m_th=original.m_th;
            m_fftb=original.m_fftb;m_fftth=original.m_fftth;
            m_dbx=original.m_dbx;m_dthx=original.m_dthx;
            m_dby=original.m_dby;m_dthy=original.m_dthy;
            m_vx=original.m_vx;
            m_vy=original.m_vy;
            m_ButterflyTraj=original.m_ButterflyTraj;
            //MAYBE CAN BE SHARED ....
            m_filter1_Lorenz63=original.m_filter1_Lorenz63;
            m_laplacienFrac_freq=original.m_laplacienFrac_freq;
            m_filter1=original.m_filter1;
            m_filter2=original.m_filter2;
            m_filter3=original.m_filter3;


            //shared constant and temporary buffers
            m_fftadv_term1.assign(original.m_fftadv_term1, true);		m_fftadv_termTh1.assign(original.m_fftadv_termTh1, true);
            m_fftadv_term2.assign(original.m_fftadv_term2, true);		m_fftadv_termTh2.assign(original.m_fftadv_termTh2, true);
            m_fftadv_term3.assign(original.m_fftadv_term3, true);		m_fftadv_termTh3.assign(original.m_fftadv_termTh3, true);
            m_fftadv_term4.assign(original.m_fftadv_term4, true);		m_fftadv_termTh4.assign(original.m_fftadv_termTh4, true);

            m_vx2.assign(original.m_vx2, true);
            m_vy2.assign(original.m_vy2, true);
            m_dbx2.assign(original.m_dbx2, true);		m_dthx2.assign(original.m_dthx2, true);
            m_dby2.assign(original.m_dby2, true);		m_dthy2.assign(original.m_dthy2, true);

            m_vx3.assign(original.m_vx3, true);
            m_vy3.assign(original.m_vy3, true);
            m_dbx3.assign(original.m_dbx3, true);		m_dthx3.assign(original.m_dthx3, true);
            m_dby3.assign(original.m_dby3, true);		m_dthy3.assign(original.m_dthy3, true);

            m_vx4.assign(original.m_vx4, true);
            m_vy4.assign(original.m_vy4, true);
            m_dbx4.assign(original.m_dbx4, true);		m_dthx4.assign(original.m_dthx4, true);
            m_dby4.assign(original.m_dby4, true);		m_dthy4.assign(original.m_dthy4, true);

            m_fx.assign(original.m_fx, true);
            m_fy.assign(original.m_fy, true);
            m_laplacien_freq.assign(original.m_laplacien_freq, true);
            m_laplacienInv_freq.assign(original.m_laplacienInv_freq, true);

            m_ddx.assign(original.m_ddx, true);
            m_ddy.assign(original.m_ddy, true);

            m_fftvTmp.assign(original.m_fftvTmp, true);
            m_fftvTmp2.assign(original.m_fftvTmp2, true);

            m_tmpL63.assign(original.m_tmpL63,true);




}

int DynamicalModel::getNx(){return m_nx;}
int DynamicalModel::getNy(){return m_ny;}
int DynamicalModel::getNt(){return m_nt;}
int DynamicalModel::getNs(){return m_nbSnapshots;}
double DynamicalModel::getTimeStep(){return m_timeStep;}
string DynamicalModel::getMainDir(){return m_mainDir;}

const CImgList<double>& DynamicalModel::getCurrentfState(){return m_fftb;}
const CImgList<double>& DynamicalModel::getCurrentfState_th(){return m_fftth;}
const CImg<double>& DynamicalModel::getCurrentState(){return m_b;}
const CImg<double>& DynamicalModel::getCurrentState_th(){return m_th;}

const CImgList<double>& DynamicalModel::getInitialfState(){return m_fftb0;}
const CImgList<double>& DynamicalModel::getInitialfState_th(){return m_fftth0;}

const CImgList<double>& DynamicalModel::getSnapshots(){return m_snap;}
int * DynamicalModel::getSnapshotIndices(){return m_idx_snap;}

void DynamicalModel::setRayleighNumber(double Ra){m_RaNumber=Ra;}

void DynamicalModel::setPrandtlNumber(double Pr){m_diffCoef=Pr;}

void DynamicalModel::setSpatialInitialCondition(const CImg<double> initialCondition){
	//b(x,0)
    m_b0=initialCondition.get_slice(0);
	m_fftb0=m_b0.get_FFT(false);
	cimglist_for(m_fftb0,l) m_fftb0[l]=m_fftb0[l]/(4.*m_nx*m_ny);
	m_zeroCross(m_fftb0);
	m_b0=m_fftb0.get_FFT(true)[0]*(4.*m_nx*m_ny);
	m_b=m_b0;
	m_fftb=m_fftb0;
	//theta(x,0)
    m_th0=initialCondition.get_slice(1);
    m_fftth0=m_th0.get_FFT(false);
	cimglist_for(m_fftth0,l) m_fftth0[l]=m_fftth0[l]/(4.*m_nx*m_ny);
	m_zeroCross(m_fftth0);
	m_th0=m_fftth0.get_FFT(true)[0]*(4.*m_nx*m_ny);
	m_th=m_th0;
	m_fftth=m_fftth0;
}

void DynamicalModel::setFreqInitialCondition(CImgList<double> initialCondition, CImgList<double> initialCondition_th){

	m_fftb0=initialCondition;
	m_zeroCross(m_fftb0);
    m_b0=m_fftb0.get_FFT(true)[0]*(4*m_nx*m_ny);
	m_b=m_b0;
	m_fftb=m_fftb0;

	m_fftth0=initialCondition_th;
	m_zeroCross(m_fftth0);
	m_th0=m_fftth0.get_FFT(true)[0]*(4*m_nx*m_ny);
	m_th=m_th0;
	m_fftth=m_fftth0;

}



const CImg<double>& DynamicalModel::setLorenzInitialCondition(double x0,double y0,double z0,double a63, double r63, double pr63){

	//sin in z on [0,M_PI] not trivial in frequency domain -> spatial initialization easier!
	m_XYZ0[0]=x0;
	m_XYZ0[1]=y0;
	m_XYZ0[2]=z0;
	m_diffCoef=pr63;
	m_a63=a63;
	m_r63=r63;
	m_RaNumber=pow((double)(m_a63*m_a63+M_PI*M_PI),3.)/(m_a63*m_a63)*m_r63;
	m_b63=4.*M_PI*M_PI/(double)(m_a63*m_a63+M_PI*M_PI);
	if ((int)(m_r63+0.5)<1){cout<<"Cautious! flow will not develop since r <1 (r,b)=("<<m_r63<<","<<m_b63<<")"<<endl;}
	else{cout<<"Flow develops with (r,b)=("<<m_r63<<","<<m_b63<<")"<<endl;}
	m_set_filters();
	return m_highDimStateFromROMLorenz63(x0,y0,z0);

}

void DynamicalModel::setAdditionalNoise(double sigma, int *frequency, int timeIndex){
	m_sigmaNoise[timeIndex]=sigma;
	m_freqNoise[timeIndex][0]=frequency[0];m_freqNoise[timeIndex][1]=frequency[1];
}
double DynamicalModel::getAdditionalNoiseStd(int timeIndex){
	return m_sigmaNoise[timeIndex];
}
int* DynamicalModel::getAdditionalNoiseFrq(int timeIndex){
	return m_freqNoise[timeIndex];
}

void DynamicalModel::display(){
	cout<<"Domain length   		= "<<m_domainLength<<endl;
	cout<<"time step       		= "<<m_timeStep<<endl;
	cout<<"Fractional exponent 	= "<<m_alpha<<endl;
	cout<<"Filter type     		= "<<m_flag_filter<<endl;
	cout<<"x-Resolution    		= "<<m_nx<<" points"<<endl;
	cout<<"y-Resolution    		= "<<m_ny<<" points"<<endl;
	cout<<"t-Resolution    		= "<<m_nt<<" points"<<endl;
	cout<<"Prandtl number  		= "<<m_diffCoef<<endl;
	cout<<"Rayleigh (r) number 		= "<<m_RaNumber<<" ("<<m_r63<<")"<<endl;
	cout<<"Saved snapshots number  = "<<m_nbSnapshots<<endl;
	if (m_a63>0){
		cout<<"Lorenz model parameter a  = "<<m_a63<<endl;
		cout<<"Lorenz model parameter b  = "<<m_b63<<endl;
	}


}
void DynamicalModel::foward(string saveDir){
	const CImgList<double> p;
	foward(saveDir,p);
}

void DynamicalModel::subspaceProj(const CImgList<double>&p, CImg<double> & b,CImg<double> & th,CImgList<double> & fftb,CImgList<double> & fftth){

	CImg<double> buffer(2.*m_nx,2.*m_ny,2,1,0.);
	buffer.get_shared_slice(0)=b;
	buffer.get_shared_slice(1)=th;
	//proj low dimensional
	CImg<double> lowDimParam(p.width(),1,1,1,0.);
	cimg_forX(lowDimParam,x){lowDimParam(x)=p[x].get_mul(buffer).sum();}
	//rec high dimensional
	buffer.fill(0.);
	for (int i=0;i<p.width();++i){buffer+=lowDimParam(i)*p[i];}
	//update spatial domain
    b=buffer.get_slice(0);
    th=buffer.get_slice(1);
	//update fourier domain
	fftb=b.get_FFT(false);fftb[0]/=(4.*m_nx*m_ny);fftb[1]/=(4.*m_nx*m_ny);
	fftth=th.get_FFT(false);fftth[0]/=(4.*m_nx*m_ny);fftth[1]/=(4.*m_nx*m_ny);

}
void DynamicalModel::foward(string saveDir,const CImgList<double>&p){


		//if (saveDir.size()==0){cout<<"Warning: not saving because path is empty ... "<<endl;}
		string saveFile;
		if (!p.is_empty()){subspaceProj(p,m_b0,m_th0,m_fftb0,m_fftth0);}
		m_b=m_b0;
		m_fftb=m_fftb0;
		m_th=m_th0;
		m_fftth=m_fftth0;
	    unsigned long clock;
	    clock = cimg::time();
		CImgList<double> DispList_b,DispList_th;
		CImgDisplay dispTmp_b(DispList_b,"vorticity ");CImgDisplay dispTmp_th(DispList_th,"temperature ");
		for (int t=0;t<=m_nt;t++){
			// save
			std::ostringstream oss; oss <<t;
			if (saveDir.size()>0){
				m_fftb.save_cimg(const_cast<char*>((m_mainDir+saveDir+"/fftb" + oss.str() +".cimg").c_str()));
				m_fftth.save_cimg(const_cast<char*>((m_mainDir+saveDir+"/fftth" + oss.str() +".cimg").c_str()));
			}
			if (t==0 || m_idx_snap[t]>0){
				m_snap[m_idx_snap[t]].get_shared_slice(0)=m_b;
				m_snap[m_idx_snap[t]].get_shared_slice(1)=m_th;
				if (saveDir.size()>0){
					m_b.save_cimg(const_cast<char*>((m_mainDir+saveDir+"/b" + oss.str() +".cimg").c_str()));
					m_th.save_cimg(const_cast<char*>((m_mainDir+saveDir+"/th" + oss.str() +".cimg").c_str()));
				}

			}
			// evaluation of the next state
			 m_foward_CrankNicolson_RungeKutta4();
			 if (!p.is_empty()){subspaceProj(p,m_b,m_th,m_fftb,m_fftth);}
			 //display
			 DispList_b=m_b;  DispList_th=m_th;
			 dispTmp_b.display(DispList_b);dispTmp_th.display(DispList_th);
			 if(t!=0){cout<<(int)((double)(t+1)/m_nt*100.-1)<<"%";fflush(stdout);}
		}
		cout << "\n one Full model foward execution time : " << (cimg::time()-clock)/1000. << " s." << endl;


}

void DynamicalModel::fowardROMLorenz63(string saveDir)
{
			//if (saveDir.size()==0){cout<<"Warning: not saving because path is empty ... "<<endl;}
			CImg<double> aux;
			for (int i=0;i<3;i++) m_XYZ[i]=m_XYZ0[i];
			aux=m_highDimStateFromROMLorenz63(m_XYZ[0],m_XYZ[1],m_XYZ[2]);
            m_b=aux.get_slice(0);m_th=aux.get_slice(1);
			//aux.display("b0/th0 lorenz63");
			CImgList<double> DispList_b,DispList_th;
			CImgDisplay dispTmp_b(DispList_b,"b Lorenz63 ");CImgDisplay dispTmp_th(DispList_th,"th Lorenz63 ");
			string saveFile;
			for (int t=0;t<=m_nt;t++){
				// save
				if (saveDir.size()>0){
					std::ostringstream oss; oss <<t;
					m_b.save_cimg(const_cast<char*>((m_mainDir+saveDir+"/b" + oss.str() +".cimg").c_str()));
					m_th.save_cimg(const_cast<char*>((m_mainDir+saveDir+"/th" + oss.str() +".cimg").c_str()));
				}
				cimg_forY(m_ButterflyTraj,l) m_ButterflyTraj(t,l)=m_XYZ[l];
				// evaluation of the next state
				m_fowardROMLorenz63_RungeKutta4();
				//display
				if (t%10==0){
					aux=m_highDimStateFromROMLorenz63(m_XYZ[0],m_XYZ[1],m_XYZ[2]);
                    m_b=aux.get_slice(0);m_th=aux.get_slice(1);
					DispList_b=m_b;DispList_th=m_th;
					dispTmp_b.display(DispList_b);dispTmp_th.display(DispList_th);
				}
				if(t!=0){cout<<(int)((double)(t+1)/m_nt*100.-1)<<"%";fflush(stdout);}

			}
			//display
			//m_ButterflyTraj.display("Butterfly Trajectory");
			// save
			(m_ButterflyTraj.get_transpose()).save_ascii(const_cast<char*>((m_mainDir+saveDir+"/butterflyLorenz63.txt").c_str()));


}

const CImgList<double>  & DynamicalModel::foward_oneStep(){
	m_foward_CrankNicolson_RungeKutta4();
	return m_fftb;
}
const CImgList<double>  & DynamicalModel::foward_oneStep(double std, int *freq, int t){
	m_foward_CrankNicolson_RungeKutta4();
	//noise generation at a given frequency
	if (freq[0]>=0 && freq[0]<m_nx && freq[1]>=0 && freq[1]<m_ny){
		m_fftb[0](freq[0],freq[1])+=std;
		m_fftb[1](freq[0],freq[1])+=std;
		if (freq[0]!=0&&freq[1]!=0){
			m_fftb[0](2*m_nx-freq[0],2*m_ny-freq[1])+=std;
			m_fftb[1](2*m_nx-freq[0],2*m_ny-freq[1])-=std;
		}
	}
	//update
	m_b=m_fftb.get_FFT(true)[0];
	m_th=m_fftth.get_FFT(true)[0];
	m_b *=4*m_nx*m_ny;
	m_th *=4*m_nx*m_ny;
	if (t==0 || m_idx_snap[t]>0){
		m_snap[m_idx_snap[t]].get_shared_slice(0)=m_b;
		m_snap[m_idx_snap[t]].get_shared_slice(1)=m_th;
	}
	if (t==1){
		m_snap[m_idx_snap[0]].get_shared_slice(0)=m_b0;
		m_snap[m_idx_snap[0]].get_shared_slice(1)=m_th0;
	}
	return m_fftb;
}


void DynamicalModel::setSaveDir(string saveDir){
	m_mainDir=saveDir;
}
string DynamicalModel::getSaveDir(){
	return m_mainDir;
}

/*Private methods*/

void DynamicalModel::m_allocate(){

	m_b0.assign(2*m_nx,2*m_ny,1,1,0.);m_th0.assign(2*m_nx,2*m_ny,1,1,0.);
	m_fftb0.assign(2,2*m_nx,2*m_ny,1,1,0.);m_fftth0.assign(2,2*m_nx,2*m_ny,1,1,0.);
	m_b.assign(2*m_nx,2*m_ny,1,1,0.);m_th.assign(2*m_nx,2*m_ny,1,1,0.);
	m_fftb.assign(2,2*m_nx,2*m_ny,1,1,0.);m_fftth.assign(2,2*m_nx,2*m_ny,1,1,0.);
	m_dbx.assign(2*m_nx,2*m_ny,1,1,0.);m_dthx.assign(2*m_nx,2*m_ny,1,1,0.);
	m_dby.assign(2*m_nx,2*m_ny,1,1,0.);m_dthy.assign(2*m_nx,2*m_ny,1,1,0.);
	m_vx.assign(2*m_nx,2*m_ny,1,1,0.);
	m_vy.assign(2*m_nx,2*m_ny,1,1,0.);

	m_filter1.assign(2*m_nx,2*m_ny,1,1,0.);
	m_filter2.assign(2*m_nx,2*m_ny,1,1,0.);
	m_filter3.assign(2*m_nx,2*m_ny,1,1,0.);

	m_fftadv_term1.assign(2,2*m_nx,2*m_ny,1,1,0.);m_fftadv_termTh1.assign(2,2*m_nx,2*m_ny,1,1,0.);
	m_fftadv_term2.assign(2,2*m_nx,2*m_ny,1,1,0.);m_fftadv_termTh2.assign(2,2*m_nx,2*m_ny,1,1,0.);
	m_fftadv_term3.assign(2,2*m_nx,2*m_ny,1,1,0.);m_fftadv_termTh3.assign(2,2*m_nx,2*m_ny,1,1,0.);
	m_fftadv_term4.assign(2,2*m_nx,2*m_ny,1,1,0.);m_fftadv_termTh4.assign(2,2*m_nx,2*m_ny,1,1,0.);

	m_vx2.assign(2*m_nx,2*m_ny,1,1,0.);
	m_vy2.assign(2*m_nx,2*m_ny,1,1,0.);
	m_dbx2.assign(2*m_nx,2*m_ny,1,1,0.);m_dthx2.assign(2*m_nx,2*m_ny,1,1,0.);
	m_dby2.assign(2*m_nx,2*m_ny,1,1,0.);m_dthy2.assign(2*m_nx,2*m_ny,1,1,0.);

	m_vx3.assign(2*m_nx,2*m_ny,1,1,0.);
	m_vy3.assign(2*m_nx,2*m_ny,1,1,0.);
	m_dbx3.assign(2*m_nx,2*m_ny,1,1,0.);m_dthx3.assign(2*m_nx,2*m_ny,1,1,0.);
	m_dby3.assign(2*m_nx,2*m_ny,1,1,0.);m_dthy3.assign(2*m_nx,2*m_ny,1,1,0.);

	m_vx4.assign(2*m_nx,2*m_ny,1,1,0.);
	m_vy4.assign(2*m_nx,2*m_ny,1,1,0.);
	m_dbx4.assign(2*m_nx,2*m_ny,1,1,0.);m_dthx4.assign(2*m_nx,2*m_ny,1,1,0.);
	m_dby4.assign(2*m_nx,2*m_ny,1,1,0.);m_dthy4.assign(2*m_nx,2*m_ny,1,1,0.);

	m_fx.assign(1,2*m_nx,1,1,0.);
	m_fy.assign(1,2*m_ny,1,1,0.);
	m_laplacien_freq.assign(2*m_nx,2*m_ny,1,1,0.);
	m_laplacienInv_freq.assign(2*m_nx,2*m_ny,1,1,0.);
	m_laplacienFrac_freq.assign(2*m_nx,2*m_ny,1,1,0.);

	m_ddx.assign(2,2*m_nx,2*m_ny,1,1,0.);
	m_ddy.assign(2,2*m_nx,2*m_ny,1,1,0.);

	m_fftvTmp.assign(2,2*m_nx,2*m_ny,1,1,0.);
	m_fftvTmp2.assign(2,2*m_nx,2*m_ny,1,1,0.);

	m_XYZ0=new double[3];
	m_XYZ=new double[3];
	m_ButterflyTraj.assign(m_nt+1,3,1,1,0.);
	m_filter1_Lorenz63.assign(3);
    m_tmpL63.assign(2.*m_nx,2.*m_ny,2,1,0.);

    m_sigmaNoise=new double[m_nt];
    m_freqNoise=new int*[m_nt];
    for (int i=0;i<m_nt;i++){
    	m_sigmaNoise[i]=0.;
    	m_freqNoise[i]=new int[2];
    	m_freqNoise[i][0]=0;m_freqNoise[i][1]=0;
    }

    m_idx_snap=new int[m_nt+1];
	m_snap.assign(m_nbSnapshots,2*m_nx,2*m_ny,2,1,0.);
}


void DynamicalModel::m_set_filters(){


	/*frequency vectors*/
    for(int k=0;k<m_nx;k++)
    {
        m_fx(0,k)=k;
    }
	for(int k=0;k<m_nx;k++){m_fx(0,k+m_nx)=-m_nx+k;}
	for(int k=0;k<m_ny;k++){m_fy(0,k)=k;}
	for(int k=0;k<m_ny;k++){m_fy(0,k+m_ny)=-m_ny+k;}

	/*laplacian filter*/
	cimg_forY(m_laplacien_freq,y){
        m_laplacien_freq.get_shared_row(y)=(m_fx.get_mul(m_fx)).get_transpose()*(-1.)*pow((double)(2.*M_PI/m_domainLength),2.);
	}
	m_laplacien_freq.transpose();

    cimg_forY(m_laplacien_freq,x)
    {
        std::cout << 'PROUT : ' << &m_fy;
        m_laplacien_freq.get_shared_row(x) += (m_fy.get_mul(m_fy)).get_transpose()*(-1.)*pow((double)(2.*M_PI/m_domainLength),2.);
	}
	m_laplacien_freq.transpose();//m_laplacien_freq.display("laplacian");

	/*minus laplacian fractional*/
	m_laplacien_freq(0,0)=-1.;
	m_laplacienFrac_freq=m_laplacien_freq.get_fill(1.).get_div((-1*m_laplacien_freq).get_pow(m_alpha));
	m_laplacienInv_freq=m_laplacien_freq.get_fill(1.).get_div((-1*m_laplacien_freq).get_pow(1.));//m_laplacienInv_freq.display("laplacian inv");
	m_laplacien_freq(0,0)=0.;

	/*derivative filters*/
	m_fx.transpose();
	cimg_forY(m_ddx[1],y){
        m_ddx[1].get_shared_row(y)=m_fx*2.*M_PI/m_domainLength;
	}
	m_fx.transpose();
	m_ddy[1].transpose();
	m_fy.transpose();
	cimg_forY(m_ddy[1],x){
            m_ddy[1].get_shared_row(x)=m_fy*2.*M_PI/m_domainLength;
	}
	m_fy.transpose();
	m_ddy[1].transpose();

	/* setup of the filters */

	if (m_flag_filter=="Shannon")
	{
		int fx_max=(int)((m_nx-1.)/2.);
		int fy_max=(int)((m_ny-1.)/2.);
		cimg_forXY(m_filter1,x,y){m_filter1(x,y)=1./(1.-0.5*m_timeStep*m_diffCoef*m_laplacien_freq(x,y));}
		cimg_forXY(m_filter2,x,y){m_filter2(x,y)=1. + 0.5*m_timeStep*m_diffCoef*m_laplacien_freq(x,y);}
		for (int x=0;x<=fx_max;x++){for (int y=0;y<=fy_max;y++){m_filter3(x,y)=1;}	}
		for (int x=2*m_nx-fx_max;x<=2*m_nx-1;x++){for  (int y=0;y<=fy_max;y++){m_filter3(x,y)=1;}}
		for (int x=0;x<=fx_max;x++){for (int y=2*m_ny-fy_max;y<=2*m_ny-1;y++){m_filter3(x,y)=1;}}
		for (int x=2*m_nx-fx_max;x<=2*m_nx-1;x++){for (int y=2*m_ny-fy_max;y<=2*m_ny-1;y++){m_filter3(x,y)=1;}	}
		double gamma=(M_PI*M_PI+m_a63*m_a63);
		m_filter1_Lorenz63(0)=1./(1./(m_timeStep)+0.5*m_diffCoef*gamma);
		m_filter1_Lorenz63(1)=1./(1./(m_timeStep)+0.5*gamma);
		m_filter1_Lorenz63(2)=m_timeStep;

	}
	else
	{cout<<"Invalid choice for m_flag_filter"<<endl;}

	m_zeroCross(m_filter1);
	m_zeroCross(m_filter2);
	m_zeroCross(m_filter3);
	m_zeroCross(m_laplacien_freq);
	m_zeroCross(m_laplacienFrac_freq);
	m_zeroCross(m_ddx[1]);
	m_zeroCross(m_ddy[1]);

}


void DynamicalModel::m_foward_CrankNicolson_RungeKutta4(){

	//evaluation advection term
	CImgList<double> fftadv_tmp(2),fftadv_tmpTh(2);
	// velocity field evaluation
	//f(c_0)
	m_velocityEval(m_fftb,m_vx,m_vy);
	m_advtermEval(m_fftb,m_fftadv_term1,m_vx,m_vy,m_dbx,m_dby);//m_fftb.display("fftb0");
	m_rayleighBernhartBTerm(m_fftth,m_fftadv_term1);
	cimglist_for(fftadv_tmp,l) fftadv_tmp[l]=m_fftb[l]+0.5*m_timeStep*m_fftadv_term1[l];
	m_advtermEval(m_fftth,m_fftadv_termTh1,m_vx,m_vy,m_dthx,m_dthy);//m_fftth.display("fftth0");m_fftadv_termTh1.display("m_fftadv_termTh1");
	m_rayleighBernhartThTerm(m_fftb,m_fftadv_termTh1);//m_fftadv_termTh1.display("m_fftadv_termTh1_");
	cimglist_for(fftadv_tmpTh,l) fftadv_tmpTh[l]=m_fftth[l]+0.5*m_timeStep*m_fftadv_termTh1[l];//m_fftadv_termTh1.display("fftadv_tmpTh");
	//f(c_1)
	m_velocityEval(fftadv_tmp,m_vx2,m_vy2);//fftadv_tmp.display("fftb1");
	m_advtermEval(fftadv_tmp,m_fftadv_term2,m_vx2,m_vy2,m_dbx2,m_dby2);
	m_rayleighBernhartBTerm(fftadv_tmpTh,m_fftadv_term2);
	cimglist_for(fftadv_tmp,l) fftadv_tmp[l]=m_fftb[l]+0.5*m_timeStep*m_fftadv_term2[l];
	m_advtermEval(fftadv_tmpTh,m_fftadv_termTh2,m_vx2,m_vy2,m_dthx2,m_dthy2);//fftadv_tmpTh.display("fftth1");
	m_rayleighBernhartThTerm(fftadv_tmp,m_fftadv_termTh2);
	cimglist_for(fftadv_tmpTh,l) fftadv_tmpTh[l]=m_fftth[l]+0.5*m_timeStep*m_fftadv_termTh2[l];
	//f(c_2)
	m_velocityEval(fftadv_tmp,m_vx3,m_vy3);
	m_advtermEval(fftadv_tmp,m_fftadv_term3,m_vx3,m_vy3,m_dbx3,m_dby3);
	m_rayleighBernhartBTerm(fftadv_tmpTh,m_fftadv_term3);
	cimglist_for(fftadv_tmp,l) fftadv_tmp[l]=m_fftb[l]+m_timeStep*m_fftadv_term3[l];
	m_advtermEval(fftadv_tmpTh,m_fftadv_termTh3,m_vx3,m_vy3,m_dthx3,m_dthy3);
	m_rayleighBernhartThTerm(fftadv_tmp,m_fftadv_termTh3);
	cimglist_for(fftadv_tmpTh,l) fftadv_tmpTh[l]=m_fftth[l]+m_timeStep*m_fftadv_termTh3[l];
	//f(c_3)
	m_velocityEval(fftadv_tmp,m_vx4,m_vy4);
	m_advtermEval(fftadv_tmp,m_fftadv_term4,m_vx4,m_vy4,m_dbx4,m_dby4);
	m_rayleighBernhartBTerm(fftadv_tmpTh,m_fftadv_term4);
	m_advtermEval(fftadv_tmpTh,m_fftadv_termTh4,m_vx4,m_vy4,m_dthx4,m_dthy4);
	m_rayleighBernhartThTerm(fftadv_tmp,m_fftadv_termTh4);

	//propagation
	cimglist_for(m_fftb,l) m_fftb[l]=m_filter1.get_mul( m_filter2.get_mul(m_fftb[l]) + m_timeStep*(1./6.)*(m_fftadv_term1[l] + 2.*m_fftadv_term2[l]+ 2.*m_fftadv_term3[l] + m_fftadv_term4[l]) );
	cimglist_for(m_fftth,l) m_fftth[l]=m_filter1.get_mul( m_filter2.get_mul(m_fftth[l]) + m_timeStep*(1./6.)*(m_fftadv_termTh1[l] + 2.*m_fftadv_termTh2[l]+ 2.*m_fftadv_termTh3[l] + m_fftadv_termTh4[l]) );
	//update
	m_b=m_fftb.get_FFT(true)[0];
	m_th=m_fftth.get_FFT(true)[0];
	m_b *=4*m_nx*m_ny;
	m_th *=4*m_nx*m_ny;
}
void  DynamicalModel::m_fowardROMLorenz63_RungeKutta4(){

		CImg<double> tmp(3),term1(3),term2(3),term3(3),term4(3);
		cimg_forX(tmp,l) tmp(l)=m_XYZ[l];
		double gamma=(M_PI*M_PI+m_a63*m_a63);
		//f(c_0)
		m_RungeKuttaEval(tmp,term1);
		cimg_forX(tmp,l) tmp(l)=m_XYZ[l]+0.5*m_timeStep*term1(l);
		//f(c_1)
		m_RungeKuttaEval(tmp,term2);
		cimg_forX(tmp,l) tmp(l)=m_XYZ[l]+0.5*m_timeStep*term2(l);
		//f(c_2)
		m_RungeKuttaEval(tmp,term3);
		cimg_forX(tmp,l) tmp(l)=m_XYZ[l]+m_timeStep*term3(l);
		//f(c_3)
		m_RungeKuttaEval(tmp,term4);
		//propagation

		cimg_forX(tmp,l) m_XYZ[l]=m_filter1_Lorenz63(l)*( m_XYZ[l]*(1./(m_timeStep)) + (1./6.)*gamma*(term1(l) + 2.*term2(l)+ 2.*term3(l) + term4(l)) );

}
const CImg<double>&  DynamicalModel::m_highDimStateFromROMLorenz63(double xx,double yy, double zz){


	double phi_t0=xx*(M_PI*M_PI+m_a63*m_a63)*sqrt(2.)/(m_a63*M_PI);
    double T1_t0=yy*sqrt(2.)/(m_r63*M_PI);
    double T2_t0=zz/(m_r63*M_PI);
	cimg_forXY(m_tmpL63,x,z){
		m_tmpL63(x,z,0)=phi_t0*sin(M_PI*1./(2*(double) m_ny)*z)*sin(m_a63*1./(2*(double) m_nx)*x);
		m_tmpL63(x,z,1)=T1_t0*sin(M_PI*1./(2*(double) m_ny)*z)*cos(m_a63*1./(2*(double) m_nx)*x)-T2_t0*sin(M_PI*1./(2*(double) m_ny)*z);
	}
	return m_ifft_fft(m_tmpL63);

}
void  DynamicalModel::m_velocityEval(const CImgList<double> &fftb_k ,CImg<double> &vx,CImg<double> &vy){

// velocity field evaluation
	m_fftvTmp.assign(fftb_k);
	cimglist_for(m_fftvTmp,l) {m_fftvTmp[l]=m_fftvTmp[l].get_mul(m_laplacienFrac_freq);}
	m_fftvTmp2[0]=   m_fftvTmp[1].get_mul(m_ddy[1]);
	m_fftvTmp2[1]=-1*m_fftvTmp[0].get_mul(m_ddy[1]);
	m_zeroCross(m_fftvTmp2);
	vx=m_fftvTmp2.get_FFT(true)[0]*4*m_nx*m_ny;

	m_fftvTmp2[0]=-1*m_fftvTmp[1].get_mul(m_ddx[1]);
	m_fftvTmp2[1]=   m_fftvTmp[0].get_mul(m_ddx[1]);
	m_zeroCross(m_fftvTmp2);
	vy=m_fftvTmp2.get_FFT(true)[0]*4*m_nx*m_ny;
}

void DynamicalModel::m_advtermEval(const CImgList<double> &fftb_k ,CImgList<double> &fftadv_term,const CImg<double> &vx,const CImg<double> &vy,CImg<double> &dbx,CImg<double> &dby){

	//  gradient of  b or th
	m_fftvTmp[0]=-1*fftb_k[1].get_mul(m_ddx[1])*4*m_nx*m_ny;
	m_fftvTmp[1]=   fftb_k[0].get_mul(m_ddx[1])*4*m_nx*m_ny;
	dbx=m_fftvTmp.get_FFT(true)[0];

	m_fftvTmp[0]=-1*fftb_k[1].get_mul(m_ddy[1])*4*m_nx*m_ny;
	m_fftvTmp[1]=   fftb_k[0].get_mul(m_ddy[1])*4*m_nx*m_ny;
	dby=m_fftvTmp.get_FFT(true)[0];

	//(vx).display("vx");(vy).display("vy");

	// advection term calculus a properly spoken
	fftadv_term=(vx.get_mul(dbx)+vy.get_mul(dby)).get_FFT(false);
	cimglist_for(fftadv_term,l) {
		fftadv_term[l].mul(m_filter3);
		fftadv_term[l]/=-4*m_nx*m_ny;
	}
	m_zeroCross(fftadv_term);
}

void DynamicalModel::m_RungeKuttaEval(const CImg<double> &tmp,CImg<double> &term){

	term(0)=(tmp(1)-tmp(0)*0.5)*m_diffCoef;
	term(1)=(m_r63*tmp(0)-tmp(0)*tmp(2)-tmp(1)*0.5);
	term(2)=(tmp(0)*tmp(1)-m_b63*tmp(2));

}
void DynamicalModel::m_rayleighBernhartBTerm(const CImgList<double> &fftth_k ,CImgList<double> &fftadv_term){

	// gradient of temperature,
	fftadv_term[0]+=-1*fftth_k[1].get_mul(m_ddx[1])*m_RaNumber*m_diffCoef;
	fftadv_term[1]+=   fftth_k[0].get_mul(m_ddx[1])*m_RaNumber*m_diffCoef;

}
void DynamicalModel::m_rayleighBernhartThTerm(const CImgList<double> &fftb_k ,CImgList<double> &fftadv_term){

	// gradient of inverse laplacian of b,
	fftadv_term[0]+=-1*(m_laplacienInv_freq.get_mul(fftb_k[1])).get_mul(m_ddx[1]);
	fftadv_term[1]+=   (m_laplacienInv_freq.get_mul(fftb_k[0])).get_mul(m_ddx[1]);

}

void DynamicalModel::m_zeroCross(CImg<double> &inOut){

	cimg_forX(inOut,x){inOut(x,m_ny)=0.;}
	cimg_forY(inOut,y){inOut(m_nx,y)=0.;}
	inOut(0,0)=0.;

}
void DynamicalModel::m_zeroCross(CImgList<double> &inOut){

	cimglist_for(inOut,l)  {m_zeroCross(inOut[l]);}
}

const CImg<double>& DynamicalModel::m_ifft_fft(const CImg<double> &in){

	CImgList<double> auxfft;
	CImg<double> auxIm;
	m_tmpL63.assign(in);
	cimg_forZ(in,z){
		auxfft.assign(in.get_shared_slice(z).get_FFT(false));
		cimglist_for(auxfft,l) auxfft[l]=auxfft[l]/(4.*m_nx*m_ny);
		m_zeroCross(auxfft);
		auxIm.assign(auxfft.get_FFT(true)[0]*(4.*m_nx*m_ny));
		cimg_forXY(m_tmpL63,x,y) m_tmpL63(x,y,z)=auxIm(x,y);
	}
	return m_tmpL63;
}




