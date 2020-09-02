#ifndef __HBT_BES_CALC_CF_CC__
#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

void CHBT_BES::CalcCF(){
	CHBT_Part *partaa,*partbb;
	int ia,ib,namax,nbmax,iq,ictheta,Nctheta=20,nsample=0;
	double qinv,psisquared,dctheta,r,ctheta;
	vector<double> x;
	x.resize(4);
	if(IDA==IDB)
		dctheta=1.0/double(Nctheta);
	else
		dctheta=2.0/double(Nctheta);
	namax=parta.size();
	nbmax=partb.size();
	printf("Beginning calculation of correlation function, namax=%d, nbmax=%d\n",namax,nbmax);
	for(ia=1;ia<namax;ia++){
		partaa=parta[ia];
		if(IDA==IDB)
			nbmax=ia;
		for(ib=0;ib<nbmax;ib++){
			if(IDA!=IDB)
				partbb=partb[ib];
			else
				partbb=parta[ib];
			qinv=Getqinv(partaa->p,partbb->p);
			if(qinv<QINVTEST){
				nsample+=Nctheta;
				for(iq=0;iq<Nqinv;iq++){
					for(ictheta=0;ictheta<Nctheta;ictheta++){
						ctheta=-1+(0.5+ictheta)*dctheta;
						CalcXR(partaa,partbb,x,r);
						psisquared=wf->CalcPsiSquared(iq,r,ctheta);
						CFqinv[iq]+=psisquared;
					}
				}
			}
		}
	}
	printf("CF for %d,%d, nsample=%d\n",IDA,IDB,nsample);
	for(iq=0;iq<Nqinv;iq++){
		CFqinv[iq]=CFqinv[iq]/double(nsample);
		printf("%6.2f %g\n",(iq+0.5)*DELqinv,CFqinv[iq]);
	}
}

void CHBT_BES::CalcCF_MC(){
	CHBT_Part *partaa,*partbb;
	int ia,ib,namax,nbmax,iq,nsample,nctheta=10,ictheta;
	double qinv,psisquared,r,ctheta;
	vector<double> x;
	x.resize(4);
	namax=parta.size();
	nbmax=partb.size();
	if(IDA==IDB)
		nbmax=namax;
	printf("Beginning calculation of correlation function, namax=%d, nbmax=%d\n",namax,nbmax);
	nsample=0;
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
			nsample+=1;
			CalcXR(partaa,partbb,x,r);
			for(iq=0;iq<Nqinv;iq++){
				for(ictheta=0;ictheta<nctheta;ictheta++){
					ctheta=-1.0+2.0*randy->ran();
					psisquared=wf->CalcPsiSquared(iq,r,ctheta);
					CFqinv[iq]+=psisquared;
				}
			}
			if((10*nsample)%NMC==0)
				printf("finished %g percent\n",100.0*nsample/double(NMC));
		}
	}
	printf("CF for %d,%d, NMC=%d\n",IDA,IDB,NMC);
	for(iq=0;iq<Nqinv;iq++){
		CFqinv[iq]=CFqinv[iq]/double(NMC*nctheta);
		printf("%6.2f %g\n",(iq+0.5)*DELqinv,CFqinv[iq]);
	}
}

/*
void CHBT_BES::CalcCF_MC_ALT(){
	CHBT_Part *partaa,*partbb;
	int ia,ib,namax,nbmax,iq,nsample;
	double qinv,psisquared,r,ctheta;
	vector<double> x;
	x.resize(4);
	namax=parta.size();
	nbmax=partb.size();
	if(IDA==IDB)
		nbmax=namax;
	printf("Beginning calculation of correlation function, namax=%d, nbmax=%d\n",namax,nbmax);
	for(iq=0;iq<Nqinv;iq++){
		printf("Calculating for iq=%d\n",iq);
		
		nsample=0;
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
			CalcXR(partaa,partbb,x,r);
			qinv=Getqinv(partaa->p,partbb->p);
			if(qinv<QINVTEST){
				ctheta=-1.0+2.0*randy->ran();
				nsample+=1;
				psisquared=wf->CalcPsiSquared(iq,r,ctheta);
				CFqinv[iq]+=psisquared;
			}
		}
	}
	printf("CF for %d,%d, NMC=%d\n",IDA,IDB,NMC);
	for(iq=0;iq<Nqinv;iq++){
		CFqinv[iq]=CFqinv[iq]/double(NMC);
		printf("%6.2f %g\n",(iq+0.5)*DELqinv,CFqinv[iq]);
	}
}
*/



#endif