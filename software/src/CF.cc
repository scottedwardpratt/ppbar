#ifndef __HBT_BES_CALC_CF_CC__
#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

int CF::NQ=50;
double CF::DELQ=2.0;
/* int CF::NPHI=4;
int CF::NRAP=3;
double CF::DELRAP=0.2;
int CF::NPT=20;
double CF::DELPT=100.0; */
CRandy* CF::randy=NULL;
CWaveFunction *CF::wf=NULL;

CF::CF(){
	cf_qinv.resize(NQ);
	cf_qout.resize(NQ);
	cf_qside.resize(NQ);
	cf_qlong.resize(NQ);
	Reset();
}
void CF::Reset(){
	nsample=0;
	for(int iq=0;iq<NQ;iq++){
		cf_qinv[iq]=cf_qout[iq]=cf_qside[iq]=cf_qlong[iq]=0.0;
	}
}

void CF::Increment(CHBT_Part *parta,CHBT_Part *partb){
	int iq,ictheta;
	const int Nctheta=6;
	double r,psisquared,ctheta;
	vector<double> x(4,0.0);
	CalcXR(parta,partb,x,r);
	if(r==r && r!=0.0){
		nsample+=1;
		for(iq=0;iq<NQ;iq++){
			for(ictheta=0;ictheta<Nctheta;ictheta++){
				ctheta=-1.0+2.0*randy->ran();
				psisquared=wf->CalcPsiSquared(iq,r,ctheta);
				if(psisquared!=psisquared){
					parta->Print();
					partb->Print();
					printf("---------- psisquared = %g ------------\n",psisquared);
					exit(1);
				}
				cf_qinv[iq]+=psisquared/double(Nctheta);
			}
			ctheta=x[1]/r;
			psisquared=wf->CalcPsiSquared(iq,r,ctheta);
			cf_qout[iq]+=psisquared;
			ctheta=x[2]/r;
			psisquared=wf->CalcPsiSquared(iq,r,ctheta);
			cf_qside[iq]+=psisquared;
			ctheta=x[3]/r;
			psisquared=wf->CalcPsiSquared(iq,r,ctheta);
			cf_qlong[iq]+=psisquared;
		}
	}
	else{
		printf("weird, r=%g\n",r);
		//parta->Print();
		//partb->Print();
		//exit(1);
	}
}

/*void CF::CalcXR(CHBT_Part *partaa,CHBT_Part *partbb,vector<double> &x,double &r){
	double dummy,gamma,gammav=0.5*((partaa->p[1]/partaa->mass)+(partbb->p[1]/partbb->mass));
	gamma=sqrt(1.0+gammav*gammav);
	for(int alpha=0;alpha<4;alpha++)
		x[alpha]=partaa->x[alpha]-partbb->x[alpha];
	dummy=x[0];
	x[0]=gamma*dummy-gammav*x[1];
	x[1]=gamma*x[1]-gammav*dummy;
	r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
}*/

void CF::CalcXR(CHBT_Part *partaa,CHBT_Part *partbb,vector<double> &x,double &r){
	vector<double> pa=(partaa->p);
	vector<double> pb=(partbb->p);
	vector<double> xa=(partaa->x);
	vector<double> xb=(partbb->x);
	vector<double> P(4);
	double M2=0.0,r2,Pdotx;
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=pa[alpha]+pb[alpha];
		x[alpha]=xa[alpha]-xb[alpha];
	}
	M2=P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3];
	Pdotx=P[0]*x[0]-P[1]*x[1]-P[2]*x[2]-P[3]*x[3];
	for(alpha=0;alpha<4;alpha++){
		x[alpha]=x[alpha]-Pdotx*P[alpha]/M2;
	}
	r2=-x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3];
	r=sqrt(r2);
	x[1]=sqrt(x[1]*x[1]-x[0]*x[0])*x[1]/fabs(x[1]);
	x[0]=0.0;
	
}

void CF::DivideByNSample(){
	int iq;
	if(nsample>0){
		double ds=double(nsample);
		for(iq=0;iq<NQ;iq++){
			cf_qinv[iq]=cf_qinv[iq]/ds;
			cf_qout[iq]=cf_qout[iq]/ds;
			cf_qside[iq]=cf_qside[iq]/ds;
			cf_qlong[iq]=cf_qlong[iq]/ds;			
		}
	}
}

void CF::Print(){
	int iq;
	printf("----- CF ------, nsample=%d\n",nsample);
	printf("q(MeV/c) CF(qinv) CF(qout) CF(side) CF(qlong)\n");
	for(iq=0;iq<NQ;iq++){
		printf("%7.3f %8.5f %8.5f %8.5f %8.5f\n",iq*DELQ,cf_qinv[iq],cf_qout[iq],cf_qside[iq],cf_qlong[iq]);
	}
}


#endif