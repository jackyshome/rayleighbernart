/*
 * main.cpp
 * ROM by particle filtering for Rayleigh-Bernart model
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#include <iostream>
#include "PffuncsRayleigh.h"

using namespace std;
using namespace cimg_library;
using namespace smc;

cv_obs * y;

int main () {

	//particle filter variable
	long lNumber = 40; //number of particles
	int dimROM=20; // dimension of the ROM's subspace
	bool pop=false;
    string respath="/Users/pheas/TemporaryBuffer/logResults.txt";

	//-launch time
	unsigned long clock;
	clock = cimg::time();
	try {
		 cout<<"\ncreating observation ..."<<endl;
		 rng *pRng=new rng[1];
		 generate_observations(pRng);


		 cout<<"\ninitializing the sampler ..."<<endl;
		 sampler<cv_state> Sampler(lNumber, SMC_HISTORY_RAM);
         moveset<cv_state> Moveset(fInitialise, fMove);
		 Sampler.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
		 Sampler.SetMoveSet(Moveset);
		 Sampler.Initialise();

		 cout<<"\nparticle filtering ..."<<endl;
		 DynamicalModel aux=Sampler.GetParticleValue(1).dmObj;
		 int nt=aux.getNt();
		 for (int t=0 ; t < nt ; ++t) {
			 Sampler.Iterate();
			 cout<<endl;
		 }
		 DynamicalModel *samples=new DynamicalModel[lNumber];
		 double *unorm_weight=new double[lNumber];
		 for (int k=0 ; k < lNumber ; ++k) {
			 samples[k]=Sampler.GetParticleValue(k).dmObj;
			 unorm_weight[k]=Sampler.GetParticleWeight(k);
		 }

		 cout<<"\nROMs computation ..."<<endl;
		 ROM Rom(samples,unorm_weight,lNumber,dimROM,0);
		 Rom.display();
		 ROM Rom_mmse(samples,unorm_weight,lNumber,dimROM,1);
		 cout<<"..."<<endl;
		 ROM Rom_bound(samples,unorm_weight,lNumber,dimROM,-1);
		 cout<<"..."<<endl;
		 if (pop) {
			 Rom.get_pop().display("pop");Rom.checkPopBasis();
			 Rom_mmse.get_pop().display("pop mmse");Rom_mmse.checkPopBasis();
			 Rom_bound.get_pop().display("pop bound");Rom_bound.checkPopBasis();
		 }
		 else{
			 Rom.get_pod().display("pod");
			 Rom_mmse.get_pod().display("pod mmse");
			 Rom_bound.get_pod().display("pod bound");
		 }
		 int *m_idx_snap=samples[0].getSnapshotIndices();
		 cout<<"snapshot indices: ";
		 for (int i=0;i<=nt;++i){cout<<m_idx_snap[i]<<" ";}
		 cout<<endl;
		 int lastSnapshotTimeIndex=nt;
		 while (m_idx_snap[lastSnapshotTimeIndex]==0){--lastSnapshotTimeIndex;}

		 cout<<"\nROM evaluation ..."<<endl;
		 CImg<double> initCond,finalCond,ApproxROM,ApproxROM_mmse,ApproxROM_bound;
		 double error, error_mmse, error_bound, meanError,meanError_mmse,meanError_bound;
		 fstream frez(respath, ios::out | ios::app);
		 if (frez){ frez  << "er"    << "\t"<< "er_mmse"    << "\t"<< "er_bound"    << "\t"  << "dimROM" << "\t" << "pop" << "\t" << "Nobs" <<  		   "\t" << "Np" << "\t"    <<      "Ns"           << "\t"  << "Nt" << "\t" << "Nx" << "\t" << "Ny" <<endl;}
		 for (int dim=dimROM;dim>0;--dim){
		   meanError=0.; meanError_mmse=0.;meanError_bound=0.;
		   for (int kObs=0;kObs<y->omObj.getNbObs();++kObs){
			 initCond.assign(y->omObj.getTruth(0,kObs));
			 if (pop) {
				 ApproxROM.assign(Rom.runPopRom(initCond, samples[0].getNs(),dim,false));
				 ApproxROM_mmse.assign(Rom_mmse.runPopRom(initCond, samples[0].getNs(),dim,false));
				 ApproxROM_bound.assign(Rom_bound.runPopRom(initCond, samples[0].getNs(),dim,false));
			 }
			 else{
				 ApproxROM.assign(Rom.runPodRom(initCond,samples[0],y->omObj.getPrandtlNumber(kObs),dim));
				 ApproxROM_mmse.assign(Rom_mmse.runPodRom(initCond,samples[0],y->omObj.getPrandtlNumber(kObs),dim));
				 ApproxROM_bound.assign(Rom_bound.runPodRom(initCond,samples[0],y->omObj.getPrandtlNumber(kObs),dim));
			 }
			 finalCond.assign(y->omObj.getTruth(lastSnapshotTimeIndex,kObs));
		 	 error=((ApproxROM-finalCond).get_abs()/(finalCond.get_abs().mean())).mean()*100;
		 	 meanError+=error;
		 	 error_mmse=((ApproxROM_mmse-finalCond).get_abs()/(finalCond.get_abs().mean())).mean()*100;
		 	 meanError_mmse+=error_mmse;
		 	 error_bound=((ApproxROM_bound-finalCond).get_abs()/(finalCond.get_abs().mean())).mean()*100;
		 	 meanError_bound+=error_bound;
		 	 cout<<" Observation No "<<kObs<<", mean abs approximation error % ="<<error<<" (using particles), "<<error_mmse<<" (using mmse estimates), "<<error_bound<<" (bound) "<<endl;
		   }
		   meanError/=y->omObj.getNbObs();
		   meanError_mmse/=y->omObj.getNbObs();
		   meanError_bound/=y->omObj.getNbObs();
		   cout<<"\n ROM dim="<<dim<<" : mean abs approximation error % (averaged over observation) ="<<meanError<<" (using particles), "<<meanError_mmse<<" (using mmse estimates), "<<meanError_bound<<" (bound) \n"<<endl;
		    if (frez){
		        frez  << error    << "\t"<< error_mmse    << "\t"<< error_bound    << "\t\t" << dim    << "\t" <<   pop << "\t" << y->omObj.getNbObs() << "\t" << lNumber << "\t" << samples[0].getNs()  << "\t"   << nt << "\t"  << initCond.width() << "\t" << initCond.height() << endl;
		    }


		 }
	}
	 catch(smc::exception  e)
	 {
		 cerr << e;
		 exit(e.lCode);
	 }
	/**/
	cout << "Total execution time : " << (cimg::time()-clock)/1000. << " s." << endl;
	cout << "****************************************************" << endl << endl;



    return 0;
}

