#ifndef __HBT_BES_CALC_CF_CC__
#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

void CHBT_BES::CalcCF_MC(){
	CHBT_Part *partaa,*partbb;
	CF *cf;
	int ia,ib,namax,nbmax;
	long long int nsample=0;
	double qinv;
	bool success=false;
	vector<double> x;
	x.resize(4);
	namax=parta.size();
	nbmax=partb.size();
	if(IDA==IDB)
		nbmax=namax;
	printf("Beginning calculation of correlation function, namax=%d, nbmax=%d, NQ=%d\n",namax,nbmax,CF::NQ);
	while(nsample<NMC){
		do{
			ia=lrint(floor(namax*randy->ran()));
			ib=lrint(floor(nbmax*randy->ran()));
		}while(ia==ib && IDA==IDB);
		partaa=parta[ia];
		if(IDA!=IDB)
			partbb=partb[ib];
		else
			partbb=parta[ib];
		qinv=Getqinv(partaa->p,partbb->p);
		if(qinv<QINVTEST){
			cf=GetCF(partaa,partbb);
			if(cf!=NULL){
				nsample+=1;
				success=true;
				cf->Increment(partaa,partbb);
			}
		}
		if((10*nsample)%NMC==0 && success==true){
			printf("finished %g percent\n",100.0*nsample/double(NMC));
			success=false;
		}
	}
	for(int irap=0;irap<NRAP;irap++){
		for(int iphi=0;iphi<NPHI;iphi++){
			for(int ipt=0;ipt<NPT;ipt++){
				CFArray[irap][iphi][ipt]->Normalize();
			}
		}
	}
}

void CHBT_BES::CalcCF_Gauss(double Rout,double Rside,double Rlong){
	int imc;
	vector<double> x;
	x.resize(4);
	x[0]=0.0;
	cfgauss->nincrement=0;
	for(imc=0;imc<NMC;imc++){
		x[1]=Rout*randy->ran_gauss();
		x[2]=Rside*randy->ran_gauss();
		x[3]=Rlong*randy->ran_gauss();
		cfgauss->Increment(x);
		if((10*(imc+1))%NMC==0){
			printf("finished %g percent\n",100.0*(imc+1)/double(NMC));
		}
	}
	cfgauss->Normalize();
}


#endif