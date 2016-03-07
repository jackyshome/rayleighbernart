/*
 * ObservationModel.cpp
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */



#include "ObservationModel.h"
ObservationModel::ObservationModel(){
    cout<<" default constructor of obs model!! "<<this<<endl;
}
ObservationModel::ObservationModel(DynamicalModel *dm,
								const int nbObs,
								int tStar,
								int tFreq,
								int SamplingRate,
								double NoiseVar,
								CImg<double> TranfertFunction):  m_nbObs(nbObs), m_tStar(tStar),m_tFreq(tFreq), m_SamplingRate(SamplingRate), m_NoiseVar(NoiseVar), m_TranfertFunction(TranfertFunction)
{
    cout<<" constructor of obs model!! "<<this<<endl;

    // checks
	if (m_nbObs<0){cout<<"Error: m_nbObs<0"<<endl;}
	if (m_tStar<0){cout<<"Error: m_tStar<0"<<endl;}
	if (m_tFreq<=0){cout<<"Error: m_tFreq<=0"<<endl;}
	if (m_SamplingRate<=0){cout<<"Error: m_SamplingRate<=0"<<endl;}
	if (m_NoiseVar<0){cout<<"Error: m_NoiseVar<0"<<endl;}

	m_dm=dm;
    m_pr=new double[m_nbObs + 1];
    m_idx_obs=new int[m_dm->getNt()+1];
    for (int i=0;i<=m_dm->getNt();i++){m_idx_obs[i]=0;}
    for (int i=tStar;i<=m_dm->getNt();i+=tFreq){m_idx_obs[i]=1;}
    cout<<"\n observation sampling index: ";
    for (int i=0;i<=m_dm->getNt();i++){cout<<m_idx_obs[i]<<" ";}
    m_obsList=new CImgList<double>[nbObs];
    for (int i=0; i<nbObs;i++){ m_obsList[i].assign(m_dm->getNt()+1);}
    m_trueList=new CImgList<double>[nbObs];
    for (int i=0; i<nbObs;i++){ m_trueList[i].assign(m_dm->getNt()+1);}
    cout<<""<<endl;
}

ObservationModel::~ObservationModel() {
    cout<<" desctructor of obs model!! "<<this<<endl;
//    delete [] m_idx_obs;
//    delete [] m_pr;
}
void ObservationModel::allocateCopy(const ObservationModel& original){
	m_nbObs=original.m_nbObs;
	m_tStar=original.m_tStar;
	m_tFreq=original.m_tFreq;
	m_SamplingRate=original.m_SamplingRate;
	m_NoiseVar=original.m_NoiseVar;
	m_TranfertFunction=original.m_TranfertFunction;
	m_dm=original.m_dm;
    m_pr=new double[m_nbObs + 1];
	for (int i=0;i<=m_nbObs;i++){m_pr[i]=original.m_pr[i];}
    m_idx_obs=new int[m_dm->getNt()+1];
	for (int i=0;i<=m_dm->getNt();i++){m_idx_obs[i]=original.m_idx_obs[i];}
	m_obsList=new CImgList<double>[m_nbObs];
	m_trueList=new CImgList<double>[m_nbObs];
	for (int i=0; i<m_nbObs;i++){
		m_obsList[i].assign(m_dm->getNt()+1);
		m_trueList[i].assign(m_dm->getNt()+1);
		for (int t=0;t<=m_dm->getNt();t++){
            m_obsList[i][t].assign(original.m_obsList[i][t]);
            m_trueList[i][t].assign(original.m_trueList[i][t]);
		}
	}
}
ObservationModel&  ObservationModel::operator= (const ObservationModel &original){
    cout<<" assignement op. of obs model!! "<<this<<"="<<&original<<endl;
    allocateCopy(original);
	return *this;
}
ObservationModel::ObservationModel(const ObservationModel& other){
    cout<<" copy constructor of obs model!! "<<this<<endl;
    allocateCopy(other);

}
void ObservationModel::generateObs(string stateDir,string saveDir, int kthObs){generateObs( stateDir, saveDir,  kthObs,-1.);}
void ObservationModel::generateObs(string stateDir,string saveDir, int kthObs, double pr){

		CImg<double> addNoise(2*m_dm->getNx(),2*m_dm->getNy(),2,1,0.),tmpobs(2*m_dm->getNx(),2*m_dm->getNy(),2,1,0.);
		CImgList<double> tmpfftb,tmpfftth;
		ostringstream oss_k;oss_k <<kthObs;
		m_pr[kthObs]=pr;
		for (int t=0;t<=m_dm->getNt();t++){
			ostringstream oss;oss <<t;
			tmpfftb.load_cimg(const_cast<char*>((m_dm->getMainDir()+stateDir+"/fftb" + oss.str()+".cimg").c_str()));
			tmpfftth.load_cimg(const_cast<char*>((m_dm->getMainDir()+stateDir+"/fftth" + oss.str()+".cimg").c_str()));
			if (m_trueList[kthObs][t].is_empty()){
				m_trueList[kthObs][t].assign(2*m_dm->getNx(),2*m_dm->getNy(),2,1);
				m_trueList[kthObs][t].get_shared_slice(0)=tmpfftb.get_FFT(true)[0]*4*m_dm->getNx()*m_dm->getNy();
				m_trueList[kthObs][t].get_shared_slice(1)=tmpfftth.get_FFT(true)[0]*4*m_dm->getNx()*m_dm->getNy();
			}
			if (m_idx_obs[t]!=0){
				cimglist_for(tmpfftb,l){
                    tmpfftb[l].mul(m_TranfertFunction.get_slice(0));
                    tmpfftth[l].mul(m_TranfertFunction.get_slice(1));
					tmpfftb[l]*=4*m_dm->getNx()*m_dm->getNy();
					tmpfftth[l]*=4*m_dm->getNx()*m_dm->getNy();
				}
				tmpfftb.FFT(true);
				tmpfftth.FFT(true);
				if (m_NoiseVar>0) addNoise.noise(m_NoiseVar,0);
                tmpobs.get_shared_slice(0)=tmpfftb[0]+addNoise.get_slice(0);
                tmpobs.get_shared_slice(1)=tmpfftth[0]+addNoise.get_slice(1);
				m_obsList[kthObs][t].assign(tmpobs);
				tmpobs.save_cimg(const_cast<char*>(((saveDir+ "/obs_"+oss_k.str()+"_step_"+oss.str()+".cimg").c_str())));
			}
		}
}
const int ObservationModel::getNbObs(){
	return m_nbObs;

}
const CImg<double> ObservationModel::getObs(int l, int kthObs){
	return m_obsList[kthObs][l];

}
const CImg<double> ObservationModel::getTruth(int l, int kthObs){
	return m_trueList[kthObs][l];

}

double ObservationModel::getPrandtlNumber(int kthObs){
	return m_pr[kthObs];
}
const CImg<double>& ObservationModel::getTF(){
	return m_TranfertFunction;

}

void ObservationModel::setTranfertFunction(const CImg<double> &TranfertFunction){

	m_TranfertFunction.assign(TranfertFunction);

}
void ObservationModel::display() {
	cout<<"First observation time = "<<m_tStar<<endl;
	cout<<"Observation frequency  = "<<m_tFreq<<endl;
	cout<<"Subsampling rate       = "<<m_SamplingRate<<endl;
	cout<<"Noise variance       = "<<m_NoiseVar<<endl;
	cout<<"time length       = "<<m_dm->getNt()<<endl;

}
