/*
 * ReducedOrderModel.cpp
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#include "PffuncsRayleigh.h"

static double precision=1.e-7;


ROM::ROM(DynamicalModel *samples, double * weight, int Np,int dim, int mmse):  m_Dim(dim),  m_mmse(mmse), m_samples(samples){

	m_Nt=samples[0].getNt();
	m_Ns=samples[0].getNs();
	m_Nx=2*samples[0].getNx();
	m_Ny=2*samples[0].getNy();

	if (mmse>=0){
		m_snaps=new CImgList<double>[Np];
		m_weight=new double[Np];
		double normCte=0;
		for (int k=0;k<Np;k++){normCte+=weight[k];}
		for (int k=0;k<Np;k++){
			m_snaps[k].assign(samples[k].getSnapshots(),false);
			m_weight[k]=weight[k]/normCte;
		}
		m_Np=Np;
	}
	if (m_mmse==-1){//replacing weighted samples by true samples.
		m_idx_snap=m_samples[0].getSnapshotIndices();
		int nbObs=y->omObj.getNbObs();
		m_snaps=new CImgList<double>[nbObs];
		m_weight=new double[nbObs];
		for (int k=0; k< nbObs;++k){
			m_snaps[k].assign(m_Ns);
			for (int t=0;t<=m_Nt;t++){
				if (m_idx_snap[t]!=0 || t==0) {m_snaps[k][m_idx_snap[t]].assign(y->omObj.getTruth(t,k));}
			}
			m_weight[k]=1./(double)nbObs;
		}
		m_Np=nbObs;
	}

	m_normalize_snapshots();

	if (m_mmse==1){//replacing weighted samples by snapshots of mmse estimates.
		CImgList<double> m_buffer_snaps(m_Ns,m_Nx,m_Ny,2,1,0.);
		for (int t=0;t<m_Ns;t++){
			for (int k=0;k<Np;k++){m_buffer_snaps[t]+=m_snaps[k][t]*m_weight[k];}
			m_snaps[0][t].assign(m_buffer_snaps[t]);
		}
		m_weight[0]=1;
		m_Np=1;
	}


}

ROM::~ROM() {
	// TODO Auto-generated destructor stub
}
void ROM::display() const{

	cout<<"\nROM parameters:"<<endl;
	(m_mmse==1)? cout<<"mmse estimates"<<endl : cout<<"number of particles="<<m_Np<<endl;
	cout<<"number of snapshots="<<m_Ns<<endl;
	cout<<"number of image columns="<<m_Nx<<endl;
	cout<<"number of image lines="<<m_Ny<<endl;
	cout<<"ROM dimension="<<m_Dim<<"\n"<<endl;

}

 void ROM::pod(){

	 if (m_Dim> m_Ns*m_Np){
		 cout<<"\nWARNING: dimension ("<<m_Dim<<") of ROM too large for POD. It is fixed to the max dimension (number of snapshots x number of particles ="<<m_Ns*m_Np<<") \n"<<endl;
		 m_Dim=m_Ns*m_Np;
	 }
	 CImgList<double> eig(m_eigen_zzt(0,m_Ns));
	 m_podBasis.assign(m_Dim);
	 for (int i=0;i<m_Dim;i++){ m_podBasis[i].assign(eig[i+1]);}
 }

CImgList<double> ROM::get_pod(){
	 pod();
	 return m_podBasis;
 }

 void ROM::pop(){

	 if (m_Dim> (m_Ns-1)*m_Np){
	 		 cout<<"\nWARNING: dimension ("<<m_Dim<<") of ROM too large for POP. It is fixed to the max dimension (number of snapshots-1 x number of particles ="<<(m_Ns-1)*m_Np<<") \n"<<endl;
	 		 m_Dim=(m_Ns-1)*m_Np;
	 	 }
 	 //--------------ccbar eigenvector computation ---------------

	 //a^*a eigendecomposition
	 CImgList<double> eig_ata(m_eigen_ztz(0,m_Ns-1));

	 //b^*b computation
	 CImg<double> btb((m_Ns-1)*m_Np,(m_Ns-1)*m_Np,1,1,0.);
	 CImg<double> Wroot(m_Np,m_Np,1,1,0.);
	 for (int k1=0;k1<m_Np;k1++){
	  		 for (int k2=0;k2<m_Np;k2++){
	  			  Wroot(k1,k2)=sqrt(m_weight[k1])*sqrt(m_weight[k2]);
	  			 for (int t1=0;t1<m_Ns-1;t1++){
	  					for (int t2=0;t2<m_Ns-1;t2++){
	  						btb(k1*(m_Ns-1)+t1,k2*(m_Ns-1)+t2)=(m_snaps[k1][t1+1].get_mul(m_snaps[k2][t2+1])).sum()*Wroot(k1,k2);
	  					}
	  			 }
	  		 }
	 }
	 //v^*b^*bv diagonalisation
	 CImg<double> vtbtbv(eig_ata[1].get_transpose()*btb*eig_ata[1]);
	 CImgList<double> Eig_vtbtbv(vtbtbv.get_symmetric_eigen());

	 //bv computation
	 CImgList<double> bv((m_Ns-1)*m_Np,m_Nx,m_Ny,2,1,0.);
	 for (int i=0;i<(m_Ns-1)*m_Np;i++){
		 bv[i].fill(0.);
	     for (int k=0;k<m_Np;k++){
	  			 for (int t=0;t<m_Ns-1;t++){
	  				bv[i]+=eig_ata[1](i,k*(m_Ns-1)+t)*m_snaps[k][t+1]*sqrt(m_weight[k]);
	  			 }
	  	  }
	 }

	 //computation of p = eigenvector of ccbar
	 CImg<double> Eig_vtbtbvXD_vtbtbv(Eig_vtbtbv[1]);
	  	 for (int j=0;j<(m_Ns-1)*m_Np;j++){
	  	 	for (int i=0;i<(m_Ns-1)*m_Np;i++){
	  	 		if (Eig_vtbtbv[0](i)>precision) {Eig_vtbtbvXD_vtbtbv(i,j)/=sqrt(Eig_vtbtbv[0](i));}
	  	 		else{Eig_vtbtbvXD_vtbtbv(i,j)=0;}
	  	 	}
	  }

	 CImgList<double> Eig_ccbar(m_Dim,m_Nx,m_Ny,2,1,0.);
	 for (int i=0;i<m_Dim;i++){
		 for (int j=0;j<(m_Ns-1)*m_Np;j++){
	  			Eig_ccbar[i]+=bv[j]*Eig_vtbtbvXD_vtbtbv(i,j);
	  	  }
	 }
	 m_popBasis.assign(Eig_ccbar);
	 if (m_Dim < (m_Ns-1)*m_Np) {m_popEigenValue.assign(Eig_vtbtbv[0].get_crop(0,0,0,m_Dim-1).get_transpose());}
	 else{m_popEigenValue.assign(Eig_vtbtbv[0].get_transpose());}

 	 //--------------cbar c eigenvector computation ---------------

	 //use of symetric eigenvalue solver only
	 CImgList<double> Eig_cbarc(m_Dim,m_Nx,m_Ny,2,1,0.);


		 CImg<double> vD((m_Ns-1)*m_Np,(m_Ns-1)*m_Np,1,1);
		 for (int i=0;i<(m_Ns-1)*m_Np;i++){
			 for (int j=0;j<(m_Ns-1)*m_Np;j++){
				 vD(i,j)=eig_ata[1](i,j)*sqrt(eig_ata[0](i));
		 	 }
		 }
		 //bvD computation
		 CImgList<double> bvD((m_Ns-1)*m_Np,m_Nx,m_Ny,2,1,0.);
		 for (int i=0;i<(m_Ns-1)*m_Np;i++){
			 bvD[i].fill(0.);
		 	 for (int k=0;k<m_Np;k++){
		 	  	for (int t=0;t<m_Ns-1;t++){
		 	  		bvD[i]+=vD(i,k*(m_Ns-1)+t)*m_snaps[k][t+1]*sqrt(m_weight[k]);
		 	  	}
		 	  }
		 }
		 CImg<double> btbvD((m_Ns-1)*m_Np,(m_Ns-1)*m_Np,1,1,0.);
		 for (int i=0;i<(m_Ns-1)*m_Np;i++){
		 	  for (int k=0;k<m_Np;k++){
		 	 		for (int t=0;t<m_Ns-1;t++){
		 	 			btbvD(i,(m_Ns-1)*k+t)=bvD[i].get_mul(m_snaps[k][t+1]*sqrt(m_weight[k])).sum();
		 	 		}
		 	 }
		 }
		 CImgList<double> abtbvD((m_Ns-1)*m_Np,m_Nx,m_Ny,2,1,0.);
		 for (int i=0;i<(m_Ns-1)*m_Np;i++){
		 	  for (int k=0;k<m_Np;k++){
		 	  	  	for (int t=0;t<m_Ns-1;t++){
		 	  	  	   abtbvD[i]+=btbvD(i,(m_Ns-1)*k+t)*(m_snaps[k][t]*sqrt(m_weight[k]));
		 	  	  	}
		 	  }
		 }
		 CImgList<double> eig_aat(m_eigen_zzt(0,m_Ns-1));
		 CImg<double> DwtabtbvD((m_Ns-1)*m_Np,(m_Ns-1)*m_Np,1,1,0.);
		 for (int i=0;i<(m_Ns-1)*m_Np;i++){
		 	 for (int j=0;j<m_Np*(m_Ns-1);j++){
		 	 	  	DwtabtbvD(i,j)=eig_aat[1+j].get_mul(abtbvD[i]).sum();
		 	 	  	if (eig_aat[0](j)>precision) {DwtabtbvD(i,j)/=eig_aat[0](j);}
		 	 	  	else{DwtabtbvD(i,j)=0.;}
		 	 }
		 }
		 //eigen decomposition of square but non symetric matrix D^-1wtabtbvD^0.5
			 CImgList<double> Eigvec_DwtabtbvD;
			 CImg<double> Eigval_DwtabtbvD;
			 m_FindEigenValAndVec_NonSym(DwtabtbvD,Eigval_DwtabtbvD,Eigvec_DwtabtbvD);
			 //Eig_vtbtbv[0].display("eigen value");Eigval_DwtabtbvD.display("idem");

			 //computation of q = eigen vector of cbarc
			 for (int j=0;j<m_Dim;j++){
			 	 	 for (int i=0;i<(m_Ns-1)*m_Np;i++){
			 	 		Eig_cbarc[j]+=eig_aat[1+i]*Eigvec_DwtabtbvD[j](i);
			 	 	 }
			 }

	 //add to popBasis and normalize q eigenvectors so that they also satisfy the coupled eigenvalue problem
	 m_popBasis.push_back(Eig_cbarc);
	 m_normalize_qBasis();

	 m_popRomBasis.assign(m_Dim,m_Dim,1,1,0.);

	 //b=q^*p computation
	 cimg_forXY(m_popRomBasis,x,y){
			 m_popRomBasis(x,y)=m_popBasis[y+m_Dim].get_mul(m_popBasis[x]).sum();
	 }
 }

CImgList<double> ROM::get_pop(){
 	 pop();
 	 return m_popBasis;
  }

CImg<double> ROM::runPopRom_HighComplexity(const CImg<double> &initCond, int Nt){

	if (m_popBasis.is_empty()){ pop();}
	m_snapsApprox.assign(Nt);
	m_snapsApprox[0].assign(initCond);
	for (int t=1;t<Nt;++t){
		m_snapsApprox[t].assign(m_fowardOneStep_pop(m_snapsApprox[t-1]));
	}
	return m_snapsApprox[Nt-1];
}
CImg<double> ROM::runPopRom(const CImg<double> &initCond, int Nt, bool saveHDTraj){
	return runPopRom(initCond,  Nt, m_Dim, saveHDTraj);
}
CImg<double> ROM::runPopRom(const CImg<double> &initCond, int Nt, int Nsubdim, bool saveHDTraj){


	if (m_popBasis.is_empty()){ pop();}


	CImgList<double> popBasisSave(m_popBasis);
	for (int i=m_Dim-Nsubdim;i>0;--i){m_popBasis.remove(m_Dim+Nsubdim+i-1);}
	for (int i=m_Dim-Nsubdim;i>0;--i){m_popBasis.remove(Nsubdim+i-1);}
	CImg<double> popRomBasisSave(m_popRomBasis);
	m_popRomBasis.assign(popRomBasisSave.get_crop(0,0,Nsubdim-1,Nsubdim-1));
	int dimSave=m_Dim;
	m_Dim=Nsubdim;


	m_snapsRomApprox.assign(Nt);
	if (saveHDTraj){m_snapsApprox.assign(Nt);}
	else{m_snapsApprox.assign(1);}

	m_snapsApprox[0].assign(initCond);
	m_snapsApprox[1].assign(m_fowardOneStep_pop(m_snapsApprox[0]));
	m_snapsRomApprox[1].assign(m_pop_projLowDim(m_snapsApprox[1]));
	for (int t=2;t<Nt;++t){
		m_snapsRomApprox[t].assign(m_fowardOneStep_popRom(m_snapsRomApprox[t-1]));
		if (saveHDTraj){m_snapsApprox[t].assign(m_pop_recHighDim(m_snapsRomApprox[t]));}
	}
	CImg<double> res(m_pop_recHighDim(m_snapsRomApprox[Nt-1]));

	m_popBasis.assign(popBasisSave);
	m_popRomBasis.assign(popRomBasisSave);
	m_Dim=dimSave;

	return res;
}
CImg<double> ROM::runPodRom(const CImg<double> &initCond,DynamicalModel &dmObj, double pr){
	return runPodRom(initCond,dmObj,pr,m_Dim);
}
CImg<double> ROM::runPodRom(const CImg<double> &initCond,DynamicalModel &dmObj, double pr, int Nsubdim){

	CImgList<double> podBasisSave(m_podBasis);
	for (int i=m_Dim-Nsubdim;i>0;--i){m_podBasis.remove(Nsubdim+i-1);}

	dmObj.setSpatialInitialCondition(initCond);
	dmObj.setPrandtlNumber(pr);
	string path;
	dmObj.foward(path,m_podBasis);
	CImg<double> res(m_Nx,m_Ny,2,1);
	res.get_shared_slice(0).assign(dmObj.getCurrentState());
	res.get_shared_slice(1).assign(dmObj.getCurrentState_th());

	m_podBasis.assign(podBasisSave);

	return res;
}


void ROM::checkPopBasis(){

	cout<<"\nPOP Rom model parameters:"<<endl;
	cout<<"\t eigenvalues "<<endl;
	cimg_forX(m_popEigenValue,x){cout<<m_popEigenValue(x)<<" ";}cout<<endl;
	cout<<"\t p connection coef"<<endl;
	for (int i=0;i<m_Dim;i++){
		 for (int j=0;j<m_Dim;j++){
			 cout<<m_popBasis[i].get_mul(m_popBasis[j]).sum()<<" ";
		 }
		 cout<<endl;
	}
	cout<<"\t q connection coef"<<endl;
	for (int i=0;i<m_Dim;i++){
		 for (int j=0;j<m_Dim;j++){
			 		 cout<<m_popBasis[i+m_Dim].get_mul(m_popBasis[j+m_Dim]).sum()<<" ";
		 }
		 cout<<endl;
	}
	cout<<"\t p,q connection coef"<<endl;
	for (int i=0;i<m_Dim;i++){
		for (int j=0;j<m_Dim;j++){
			cout<<m_popBasis[i+m_Dim].get_mul(m_popBasis[j]).sum()<<" ";
		}
		cout<<endl;
	}
	cout<<"\nChecks:"<<endl;

	cout<<"\t p_j^*q_j (diagonal of matrix of connection coef )"<<endl;

	CImg<double> btqj2(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
	for (int i=0;i<m_Dim;i++){
    for (int k=0;k<m_Np;k++){
 			 for (int t=0;t<m_Ns-1;t++){
 				btqj2(i,(m_Ns-1)*k+t)=m_popBasis[i+m_Dim].get_mul(m_snaps[k][t+1]*sqrt(m_weight[k])).sum();
 			 }
 	  }
	}
	CImg<double> qj1ta((m_Ns-1)*m_Np,m_Dim,1,1,0.);
	for (int i=0;i<m_Dim;i++){
    for (int k=0;k<m_Np;k++){
 			 for (int t=0;t<m_Ns-1;t++){
 				qj1ta((m_Ns-1)*k+t,i)=m_popBasis[i+m_Dim].get_mul(m_snaps[k][t]*sqrt(m_weight[k])).sum()/m_popEigenValue(i);
 			 }
 	  }
	}
    CImg<double> qj1tCtqj2(qj1ta*btqj2);
    for (int i=0;i<m_Dim;i++){
    		//for (int j=0;j<m_Dim;j++){
    			cout<<qj1tCtqj2(i,i)<<" ";
    		//}
    		//	cout<<endl;
    }
    cout<<endl;
	CImg<double> btpj2(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
		for (int i=0;i<m_Dim;i++){
	    for (int k=0;k<m_Np;k++){
	 			 for (int t=0;t<m_Ns-1;t++){
	 				btpj2(i,(m_Ns-1)*k+t)=m_popBasis[i].get_mul(m_snaps[k][t+1]*sqrt(m_weight[k])).sum();
	 			 }
	 	  }
		}
		cout<<"\t qjtctpj/lmbda  (should be a one vector)"<<endl;
	 CImg<double> qj1tCtpj2(qj1ta*btpj2);
	 for (int i=0;i<m_Dim;i++){cout<<qj1tCtpj2(i,i)<<" ";}
	  cout<<endl;

	  cout<<"\t check cbar pj -qj =0 (mean abs error %):"<<endl;
	  CImgList<double> abtp(m_Dim,m_Nx,m_Ny,2,1,0.);
	  for (int i=0;i<m_Dim;i++){
	  	 for (int k=0;k<m_Np;k++){
	  	  	 	for (int t=0;t<m_Ns-1;t++){
	  	  	 			abtp[i]+=btpj2(i,(m_Ns-1)*k+t)*(m_snaps[k][t]*sqrt(m_weight[k]));
	  	  	 	 }
	  	  }
	  }
	  CImgList<double> eig_aat(m_eigen_zzt(0,m_Ns-1));
	  CImg<double> Dwtabtp(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
	  for (int i=0;i<m_Dim;i++){
	 	  	 for (int j=0;j<m_Np*(m_Ns-1);j++){
	 	  		Dwtabtp(i,j)=eig_aat[1+j].get_mul(abtp[i]).sum();
	 	  		if (eig_aat[0](j)>precision) {Dwtabtp(i,j)/=eig_aat[0](j);}
	 	  		else{Dwtabtp(i,j)=0.;}
	 	  	 }
	  }
	  CImgList<double> wDwtabtp(m_Dim,m_Nx,m_Ny,2,1,0.);
	  for (int i=0;i<m_Dim;i++){
		  for (int j=0;j<m_Np*(m_Ns-1);j++){
			  wDwtabtp[i]+=eig_aat[1+j]*Dwtabtp(i,j);
		  }
	  }
	  CImgList<double> disp(m_Dim,m_Nx,m_Ny,2,1,0.);
	  for (int i=0;i<m_Dim;i++){disp[i]=(wDwtabtp[i]-m_popBasis[i+m_Dim]).get_abs()/(m_popBasis[i+m_Dim].get_abs().mean());}
	 // wDwtabtp.display("q=cbarp");
	  //m_popBasis.get_remove(0,disp.size()-1).display("computed q");
	  //disp.display("relative abs error |cpj-qj|/|qj|");
	  for (int i=0;i<m_Dim;i++){cout<<disp[i].mean()*100<<", ";}cout<<endl;

	  cout<<"\t check ccbar pj -lmbda pj =0 (mean abs error %):"<<endl;
	  CImg<double> atwDwtabtp(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
	  for (int i=0;i<m_Dim;i++){
	  	 for (int k=0;k<m_Np;k++){
	  	  	 	for (int t=0;t<m_Ns-1;t++){
	  	  	 		atwDwtabtp(i,(m_Ns-1)*k+t)=wDwtabtp[i].get_mul(m_snaps[k][t]*sqrt(m_weight[k])).sum();
	  	  	 	 }
	  	  }
	  }
	  CImgList<double> batwDwtabtp(m_Dim,m_Nx,m_Ny,2,1,0.);
	  for (int i=0;i<m_Dim;i++){
	  	  	 for (int k=0;k<m_Np;k++){
	  	  	  	 	for (int t=0;t<m_Ns-1;t++){
	  	  	  	 	batwDwtabtp[i]+=(m_snaps[k][t+1]*sqrt(m_weight[k]))*atwDwtabtp(i,(m_Ns-1)*k+t);
	  	  	  	 	 }
	  	  	  }
	  }
	  for (int i=0;i<m_Dim;i++){batwDwtabtp[i]/=m_popEigenValue[i];}
	  //batwDwtabtp.display("p=ccbarp/lmbda");
	  //m_popBasis.get_remove(batwDwtabtp.size(),m_popBasis.size()-1).display("computed p");
	  for (int i=0;i<m_Dim;i++){batwDwtabtp[i]=(batwDwtabtp[i]-m_popBasis[i]).get_abs()/(m_popBasis[i].get_abs().mean());}
	  //batwDwtabtp.display("relative abs error |ccbarpj-pj|/|pj|");
	  for (int i=0;i<m_Dim;i++){cout<<batwDwtabtp[i].mean()*100<<", ";}cout<<endl;


	  cout<<"\t check c qj -lmbda pj =0 (mean abs error %):"<<endl;

	  	  CImg<double> atq(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
	  	  for (int i=0;i<m_Dim;i++){
	  	  	 for (int k=0;k<m_Np;k++){
	  	  	  	 	for (int t=0;t<m_Ns-1;t++){
	  	  	  	 		atq(i,(m_Ns-1)*k+t)=m_popBasis[i+m_Dim].get_mul(m_snaps[k][t]*sqrt(m_weight[k])).sum();
	  	  	  	 	 }
	  	  	  }
	  	  }
	  	  CImgList<double> batq(m_Dim,m_Nx,m_Ny,2,1,0.);
	  	  for (int i=0;i<m_Dim;i++){
	  	  	  	 for (int k=0;k<m_Np;k++){
	  	  	  	  	 for (int t=0;t<m_Ns-1;t++){
	  	  	  	  	 	batq[i]+=atq(i,(m_Ns-1)*k+t)*(m_snaps[k][t+1]*sqrt(m_weight[k]));
	  	  	  	  	 }
	  	  	  	  }
	  	  }
	  	for (int i=0;i<m_Dim;i++){disp[i]=batq[i]/m_popEigenValue[i];}
	  	//disp.display("p=cq/lmbda");
	    //m_popBasis.get_remove(disp.size(),m_popBasis.size()-1).display("computed p");
	   for (int i=0;i<m_Dim;i++){disp[i]=(disp[i]-m_popBasis[i]).get_abs()/(m_popBasis[i].get_abs().mean());}
	   //disp.display("relative abs error |cqj/lmbda-pj|/|pj| ");
	   for (int i=0;i<m_Dim;i++){cout<<disp[i].mean()*100<<", ";}cout<<endl;


	  	cout<<"\t check cbarc qj -lmbda qj =0 (mean abs error %):"<<endl;

	  	CImg<double> btbatq(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
	  	 for (int i=0;i<m_Dim;i++){
	  		  	  for (int k=0;k<m_Np;k++){
	  		  	  	  	 	for (int t=0;t<m_Ns-1;t++){
	  		  	  	  	 		btbatq(i,(m_Ns-1)*k+t)=batq[i].get_mul(m_snaps[k][t+1]*sqrt(m_weight[k])).sum();
	  		  	  	  	 	 }
	  		  	  }
	  	 }
	  	CImgList<double> abtbatq(m_Dim,m_Nx,m_Ny,2,1,0.);
	  		  for (int i=0;i<m_Dim;i++){
	  		  	 for (int k=0;k<m_Np;k++){
	  		  	  	 	for (int t=0;t<m_Ns-1;t++){
	  		  	  	 	abtbatq[i]+=btbatq(i,(m_Ns-1)*k+t)*(m_snaps[k][t]*sqrt(m_weight[k]));
	  		  	  	 	 }
	  		  	  }
	  		  }
	  		  CImg<double> Dwtabtbatq(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
	  		  for (int i=0;i<m_Dim;i++){
	  		 	  	 for (int j=0;j<m_Np*(m_Ns-1);j++){
	  		 	  	Dwtabtbatq(i,j)=eig_aat[1+j].get_mul(abtbatq[i]).sum();
	  		 	  		if (eig_aat[0](j)>precision) {Dwtabtbatq(i,j)/=eig_aat[0](j);}
	  		 	  		else{Dwtabtbatq(i,j)=0.;}
	  		 	  	 }
	  		  }
	  		  CImgList<double> wDwtabtbatq(m_Dim,m_Nx,m_Ny,2,1,0.);
	  		  for (int i=0;i<m_Dim;i++){
	  			  for (int j=0;j<m_Np*(m_Ns-1);j++){
	  				wDwtabtbatq[i]+=Dwtabtbatq(i,j)*eig_aat[1+j];
	  			  }
	  		  }
	  		  for (int i=0;i<m_Dim;i++){wDwtabtbatq[i]/=m_popEigenValue[i];}
	  		  //wDwtabtbatq.display("q=cbarcq/lmbda");
	  		  //m_popBasis.get_remove(0,wDwtabtbatq.size()-1).display("computed q");
	  		  for (int i=0;i<m_Dim;i++){wDwtabtbatq[i]=(wDwtabtbatq[i]-m_popBasis[i+m_Dim]).get_abs()/(m_popBasis[i+m_Dim].get_abs().mean());}
	  		  //batwDwtabtp.display("relative abs error |cbarcqj-qj|/|qj|");
	  		 for (int i=0;i<m_Dim;i++){cout<<wDwtabtbatq[i].mean()*100<<", ";}cout<<endl;

}

//private methods
void ROM::m_normalize_qBasis(){



	    cout<<"Normalizing  qj's s.t. (cbar pj - qj) =0 ..."<<endl;
		CImg<double> btp(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
			for (int i=0;i<m_Dim;i++){
		    for (int k=0;k<m_Np;k++){
		 			 for (int t=0;t<m_Ns-1;t++){
		 				btp(i,(m_Ns-1)*k+t)=m_popBasis[i].get_mul(m_snaps[k][t+1]*sqrt(m_weight[k])).sum();
		 			 }
		 	  }
			}
			CImgList<double> abtp(m_Dim,m_Nx,m_Ny,2,1,0.);
		  for (int i=0;i<m_Dim;i++){
		  	 for (int k=0;k<m_Np;k++){
		  	  	 	for (int t=0;t<m_Ns-1;t++){
		  	  	 			abtp[i]+=btp(i,(m_Ns-1)*k+t)*(m_snaps[k][t]*sqrt(m_weight[k]));
		  	  	 	 }
		  	  }
		  }
		  CImgList<double> eig_aat(m_eigen_zzt(0,m_Ns-1));
		  CImg<double> Dwtabtp(m_Dim,(m_Ns-1)*m_Np,1,1,0.);
		  for (int i=0;i<m_Dim;i++){
		 	  	 for (int j=0;j<m_Np*(m_Ns-1);j++){
		 	  		Dwtabtp(i,j)=eig_aat[1+j].get_mul(abtp[i]).sum();
		 	  		if (eig_aat[0](j)>precision) {Dwtabtp(i,j)/=eig_aat[0](j);}
		 	  		else{Dwtabtp(i,j)=0.;}
		 	  	 }
		  }
		  CImgList<double> wDwtabtp(m_Dim,m_Nx,m_Ny,2,1,0.);
		  for (int i=0;i<m_Dim;i++){
			  for (int j=0;j<m_Np*(m_Ns-1);j++){
				  wDwtabtp[i]+=eig_aat[1+j]*Dwtabtp(i,j);
			  }
		  }
		  double gamma;
		  for (int i=0;i<m_Dim;i++){
			  gamma=wDwtabtp[i].get_div(m_popBasis[i+m_Dim]).mean();
			  m_popBasis[i+m_Dim]*=gamma;
		  }


}

void ROM::m_normalize_snapshots(){

	double	normalizing_cst=0;
	int count=0;
	CImg<double> W(m_Np,1,1,1,0.);
 	for (int t1=0;t1<m_Ns;t1++){
 				 for (int k1=0;k1<m_Np;k1++){
		 		 	 W(k1)=(m_weight[k1])*(m_weight[k1]);
		 			 normalizing_cst+=(m_snaps[k1][t1].get_mul(m_snaps[k1][t1])).sum()*W(k1);
		 		 }
 				 ++count;

	}
	normalizing_cst=sqrt(normalizing_cst/count);
	for (int t1=0;t1<m_Ns;t1++){
		 	 		for (int k1=0;k1<m_Np;k1++){
		 	 			m_snaps[k1][t1]/=normalizing_cst;
		 	 		}
	}
}


 void ROM::m_FindEigenValAndVec_NonSym(const CImg<double> &A, CImg<double> &valEigen,CImgList<double> &vecEigen_){

 	int size=A.width();
	valEigen.assign(size,1,1,1,0.);
	vecEigen_.assign(size,size,1,1,1,0.);
	if (size<=2){
 		CImg <double> vecEigen(size,size);
 		A.eigen(valEigen,vecEigen);
 		for (int k=0;k<size;k++){
 			for (int l=0;l<size;l++){
 				vecEigen_[k](l)=vecEigen(k,l);
 			}
 		}

 	}
 	else{
 		char jobvl='N';
 		char jobvr='V';
 		int  /*LWORK = size*size,*/ INFO;
 		double * AA=new double [size*size];
 		double * VL=new double [size*size];
 		double * VR=new double [size*size];
 		double * WR=new double [size];
 		double * WI =new double [size];
 		double *WORK;// = new double[LWORK];
 		cimg_forXY(A,k,l) AA[k*size+l] = A(k,l);
 		/* Query and allocate the optimal workspace */
 		int    LWORK = -1;
 		double wkopt;
 		dgeev_(&jobvl,&jobvr,&size, AA ,&size, WR,WI, VL, &size, VR,&size,&wkopt, &LWORK, &INFO );
 		LWORK = (int)wkopt;
 		WORK = (double*)malloc( LWORK*sizeof(double) );
 	    /* Solve eigenproblem */
 		dgeev_(&jobvl,&jobvr,&size, AA ,&size, WR,WI, VL, &size, VR,&size,WORK,&LWORK, &INFO);
 		if (INFO)
 			cout<<"eigen() : LAPACK library function dgeev_() returned error code "<< INFO<<endl;
         if (!INFO) {
 			for  (int i=0;i<size;i++){	valEigen(i) = WR[size-1-i];
 				for  (int k=0;k<size;k++){	 vecEigen_[k](i)= (VR[(size-1-k)*size+i]);}
 			}
         }
         vecEigen_.reverse();
         valEigen.mirror('x');

 		delete [] AA;delete [] VL;delete [] VR;delete [] WR;delete [] WI;delete [] WORK;
 	}
 }

 CImgList<double> ROM::m_eigen_ztz( int timeIndexMin,int Ns){

	     //x^*x computation
		 CImg<double> xtx(Ns*m_Np,Ns*m_Np,1,1,0.);
	 	 CImg<double> Wroot(m_Np,m_Np,1,1,0.);
	 	 for (int k1=0;k1<m_Np;k1++){
	 		 for (int k2=0;k2<m_Np;k2++){
	 			  Wroot(k1,k2)=sqrt(m_weight[k1])*sqrt(m_weight[k2]);
	 			 for (int t1=0;t1<Ns;t1++){
	 					for (int t2=0;t2<Ns;t2++){
	 						xtx(k1*Ns+t1,k2*Ns+t2)=(m_snaps[k1][t1+timeIndexMin].get_mul(m_snaps[k2][t2+timeIndexMin])).sum()*Wroot(k1,k2);
	 					}
	 			 }
	 		 }
	 	 }

	 	 //x^*x diagonalization
	 	return xtx.get_symmetric_eigen();
 }

 CImgList<double> ROM::m_eigen_zzt_from_ztz(const CImgList<double>& eigen_ztz,  int timeIndexMin,int Ns){

	    //xx^* eigenvector computation form x^*x eigendecomposition
	  	 CImg<double> vXdiag(eigen_ztz[1]);
	  	 for (int j=0;j<m_Np*Ns;j++){
	  	 	for (int i=0;i<m_Np*Ns;i++){
	  	 		if (eigen_ztz[0](i)>precision) {vXdiag(i,j)/=sqrt(eigen_ztz[0](i));}
	  	 		else{vXdiag(i,j)=0;}
	  	 	}
	  	 }
	  	CImgList<double>  m_aux(m_Np*Ns+1);
	  	m_aux[0].assign(eigen_ztz[0]);
	  	CImg<double> buffer(m_Nx,m_Ny,2,1,0.);
	  	for (int i=0;i<m_Np*Ns;i++){
	  		 buffer.fill(0.);
	  		 for (int k=0;k<m_Np;k++){
	  			 for (int t=0;t<Ns;t++){
	  				 buffer+=vXdiag(i,k*Ns+t)*m_snaps[k][t+timeIndexMin]*sqrt(m_weight[k]);
	  			 }
	  		 }
	  		m_aux[i+1].assign(buffer);
	  	 }
	  	 return  m_aux;

 }

 CImgList<double> ROM::m_eigen_zzt(int timeIndexMin,int Ns){

 	 //x^*x computation
	 CImgList<double> EigDec(m_eigen_ztz( timeIndexMin, Ns));

 	 //xx^* eigenvector computation
	 return m_eigen_zzt_from_ztz(EigDec, timeIndexMin, Ns);

  }
 CImg<double> ROM::m_fowardOneStep_pop(const CImg<double> &buffer){

	 CImg<double> res(buffer);
	 CImg<double> aux(m_Dim);
	 for (int i=0; i<  m_Dim;++i){
		 aux(i)=m_popBasis[i+m_Dim].get_mul(buffer).sum();
	 }
	 res.fill(0.);
	 for (int i=0; i<  m_Dim;++i){
		 res+=aux(i)*m_popBasis[i];
	 }
  	 return res;
  }

 CImg<double> ROM::m_fowardOneStep_popRom(const CImg<double> &buffer){
	 return (m_popRomBasis*(buffer.get_transpose())).get_transpose();
 }
 CImg<double> ROM::m_pop_projLowDim(const CImg<double> &HighDimState){

 	 CImg<double> res(m_Dim,1,1,1,0.);
	 cimg_forX(res,x){
		 res(x)=m_popBasis[x].get_mul(HighDimState).sum();
	 }
 	 return res;

  }
 CImg<double> ROM::m_pop_recHighDim(const CImg<double> &lowDimParam){

	 CImg<double> res(m_popBasis[0],"xyzc",0.);
	 for (int i=0;i<m_Dim;++i){
		 res+=lowDimParam(i)*m_popBasis[i];
	 }
 	 return res;


  }
