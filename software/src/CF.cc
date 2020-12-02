#ifndef __HBT_BES_CALC_CF_CC__
#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

int CF::NQ=40;
double CF::DELQ=2.0;
CHBT_BES *CF::hbt=NULL;
double CF::OUTSIDELONG_DIRECTION_CUT=0.9;
double CF::OUTSIDELONG_Q_CUT=10.0;
bool CF::USE_OUTSIDELONG_Q_CUT=true;
bool CF::USE_OUTSIDELONG_DIRECTION_CUT=false;
double CF::Dxyz=0.5;
bool CF::COAL_USE_WF=false;
vector<double> CF::psi_coal={};
int CF::Nxyz=100;
double CF::Rcoalescence=1.5;
double CF::COAL_DELR=0.05;
CRandy* CF::randy=NULL;
int CF::NSAMPLE_THETAPHI=4;

CF::CF(){
	cf_qinv.resize(NQ);
	cf_qout.resize(NQ);
	cf_qside.resize(NQ);
	cf_qlong.resize(NQ);
	norm_qinv.resize(NQ);
	norm_qout.resize(NQ);
	norm_qside.resize(NQ);
	norm_qlong.resize(NQ);
	D3q=1.0;
	Reset();
}

void CF::Reset(){
	nincrement=ncoalescence=0;
	for(int iq=0;iq<NQ;iq++){
		cf_qinv[iq]=cf_qout[iq]=cf_qside[iq]=cf_qlong[iq]=0.0;
		norm_qinv[iq]=norm_qout[iq]=norm_qside[iq]=norm_qlong[iq]=0.0;
	}
}

void CF::Increment(CHBT_Part *parta,CHBT_Part *partb){
	int iq,ithetaphi;
	double phi,ctheta,stheta,qx,qy,qz,qinv;
	double r,psisquared,ctheta_qr;
	vector<double> x(4,0.0);
	CalcXR(parta,partb,x,r);
	//x[1]=7.0*randy->ran_gauss();
	//x[2]=5.0*randy->ran_gauss();
	//x[3]=5.0*randy->ran_gauss();
	//r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
	nincrement+=1;
	if(r==r && r!=0.0){
		ncoalescence+=CoalescenceWeight(r);
		
		// Increment correlation function
		for(ithetaphi=0;ithetaphi<NSAMPLE_THETAPHI;ithetaphi++){
			phi=2.0*PI*randy->ran();
			ctheta=-1.0+2.0*randy->ran();
			stheta=sqrt(1.0-ctheta*ctheta);
			qz=ctheta;
			qx=stheta*cos(phi);
			qy=stheta*sin(phi);
			ctheta_qr=(qx*x[1]+qy*x[2]+qz*x[3])/r;
			if(fabs(ctheta_qr)>1.0){
				printf("ctheta_qr out of range\n");
				exit(1);
			}		
			for(iq=0;iq<NQ;iq++){
				qinv=(iq+0.5)*DELQ;
				psisquared=wf->CalcPsiSquared(iq,r,ctheta_qr);
				//psisquared=1.0;
				if(psisquared!=psisquared){
					parta->Print();
					partb->Print();
					printf("---------- psisquared = %g ------------\n",psisquared);
					exit(1);
				}
				cf_qinv[iq]+=psisquared/D3q;
				norm_qinv[iq]+=1.0/D3q;
				/*
				if(USE_OUTSIDELONG_DIRECTION_CUT){
				if(fabs(qx)>OUTSIDELONG_DIRECTION_CUT){
				cf_qout[iq]+=psisquared;
				norm_qout[iq]+=1;
				}
				if(fabs(qy)>OUTSIDELONG_DIRECTION_CUT){
				cf_qside[iq]+=psisquared;
				norm_qside[iq]+=1;
				}
				if(fabs(qz)>OUTSIDELONG_DIRECTION_CUT){
				cf_qlong[iq]+=psisquared;
				norm_qlong[iq]+=1;
				}
				}*/
				if(USE_OUTSIDELONG_Q_CUT){
					double qperp;
					qperp=qinv*sqrt(1.0-qx*qx);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qout[iq]+=psisquared/D3q;
						norm_qout[iq]+=1.0/D3q;
					}
					qperp=qinv*sqrt(1.0-qy*qy);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qside[iq]+=psisquared/D3q;
						norm_qside[iq]+=1.0/D3q;
					}
					qperp=qinv*sqrt(1.0-qz*qz);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qlong[iq]+=psisquared/D3q;
						norm_qlong[iq]+=1.0/D3q;
					}
				}
			}
		}
	}
	else{
		printf("weird, r=%g\n",r);
		exit(1);
		//parta->Print();
		//partb->Print();
		//exit(1);
	}
}

void CF::Increment(vector<double> &x){

	int iq,ithetaphi;
	double phi,ctheta,stheta,qx,qy,qz,qinv,qperp;
	double r,psisquared,ctheta_qr;
	r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
	nincrement+=1;

	for(ithetaphi=0;ithetaphi<NSAMPLE_THETAPHI;ithetaphi++){
		phi=2.0*PI*randy->ran();
		ctheta=-1.0+2.0*randy->ran();
		stheta=sqrt(1.0-ctheta*ctheta);
		qz=ctheta;
		qx=stheta*cos(phi);
		qy=stheta*sin(phi);
		ctheta_qr=(qx*x[1]+qy*x[2]+qz*x[3])/r;
		if(fabs(ctheta_qr)>1.0){
			printf("ctheta_qr out of range\n");
			exit(1);
		}		
		for(iq=0;iq<NQ;iq++){
			qinv=(iq+0.5)*DELQ;
			psisquared=wf->CalcPsiSquared(iq,r,ctheta_qr);
			//psisquared=1.0;
			cf_qinv[iq]+=psisquared/D3q;
			norm_qinv[iq]+=1.0/D3q;
			if(USE_OUTSIDELONG_Q_CUT){
				qperp=qinv*sqrt(1.0-qx*qx);
				if(qperp<OUTSIDELONG_Q_CUT){
					cf_qout[iq]+=psisquared/D3q;
					norm_qout[iq]+=1.0/D3q;
				}
				qperp=qinv*sqrt(1.0-qy*qy);
				if(qperp<OUTSIDELONG_Q_CUT){
					cf_qside[iq]+=psisquared/D3q;
					norm_qside[iq]+=1.0/D3q;
				}
				qperp=qinv*sqrt(1.0-qz*qz);
				if(qperp<OUTSIDELONG_Q_CUT){
					cf_qlong[iq]+=psisquared/D3q;
					norm_qlong[iq]+=1.0/D3q;
				}
			}
		}
	}

}

void CF::CalcXR(CHBT_Part *partaa,CHBT_Part *partbb,vector<double> &x,double &r){
	vector<double> P(4);
	double M2=0.0,r2,Pdotx;
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=partaa->p[alpha]+partbb->p[alpha];
		x[alpha]=partaa->x[alpha]-partbb->x[alpha];
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

void CF::Normalize(){
	int iq;
	for(iq=0;iq<NQ;iq++){
		if(norm_qinv[iq]>0){
			cf_qinv[iq]=cf_qinv[iq]/norm_qinv[iq];
		}
		else
			cf_qinv[iq]=0;
		if(norm_qout[iq]>1.0E-20){
			cf_qout[iq]=cf_qout[iq]/norm_qout[iq];
		}
		else
			cf_qout[iq]=0;
		if(norm_qside[iq]>1.0E-20){
			cf_qside[iq]=cf_qside[iq]/norm_qside[iq];
		}
		else
			cf_qside[iq]=0;
		if(norm_qlong[iq]>1.0E-20){
			cf_qlong[iq]=cf_qlong[iq]/norm_qlong[iq];
		}
		else
			cf_qlong[iq]=0;
	}
}

void CF::Print(){
	int iq;
	printf("#---- CF ------, norm_qinv=%g\n",norm_qinv[0]);
	printf("q(MeV/c) CF(qinv) CF(qout) CF(side) CF(qlong)\n");
	for(iq=0;iq<NQ;iq++){
		printf("%7.3f %8.5f %8.5f %8.5f %8.5f\n",(iq+0.5)*DELQ,cf_qinv[iq],cf_qout[iq],cf_qside[iq],cf_qlong[iq]);
	}
	printf("nincrement=%lld\n",nincrement);
}

double CF::CoalescenceWeight(double r){
	double weight;
	int nr=psi_coal.size();
	if(COAL_USE_WF){
		int ir;
		ir=floorl(r/COAL_DELR);
		if(ir<nr){
			weight=psi_coal[ir]*psi_coal[ir]/D3q;
		}
		else
			weight=0.0;
	}
	else{
		weight=exp(-0.5*r*r/(Rcoalescence*Rcoalescence));
		weight=weight/(D3q*pow(2.0*PI*Rcoalescence*Rcoalescence,1.5));
	}
	return weight;
}

void CF::CalcCoalWF(){
	COAL_DELR=0.05;
	const int nr=1000;
	double root4pi=sqrt(4.0*PI);
	int ir,ntries=0;
	psi_coal.resize(nr);
	double A,C;   // prefactors
	double B=2.24;   // binding energy
	double mu=0.5*939.0;  // reduced mass
	double a=1.4;    // well width
	double k,q;  // wave number inside/outside well
	double f,dfdk,dk,norm,r;
	q=sqrt(2.0*mu*B)/HBARC;
	k=2.0/a;   // a guess
	f=1.0+(q/k)*tan(k*a);
	do{
		ntries+=1;
		dfdk=-(q/(k*k))*tan(k*a)+(q*a/k)/pow(cos(k*a),2);
		dk=-f/dfdk;
		if(fabs(dk)>0.05)
			dk=0.05*dk/fabs(dk);
		k+=dk;
		f=1.0+(q/k)*tan(k*a);
		//printf("k=%g, f=%g\n",k,f);
	}while(fabs(f)>1.0E-12 && ntries<100);
	if(ntries==100){
		printf("Coalescence WF calculation failed to converge, will exit\n");
		exit(1);
	}
	A=1.0/(sin(k*a));
	norm=A*A*((0.5*a)-(0.25/k)*sin(2.0*k*a));
	norm+=0.5/q;
	A=A/sqrt(norm);
	C=1.0/sqrt(norm);
	for(ir=0;ir<nr;ir++){
		r=(ir+0.5)*COAL_DELR;
		if(r<a)
			psi_coal[ir]=A*sin(k*r);
		else
			psi_coal[ir]=C*exp(-q*(r-a));
		psi_coal[ir]=psi_coal[ir]/(root4pi*r);
	}
	double normcheck=0.0,r2bar=0.0;
	for(ir=0;ir<nr;ir++){
		r=(0.5+ir)*COAL_DELR;
		r2bar+=COAL_DELR*4.0*PI*r*r*psi_coal[ir]*psi_coal[ir]*r*r;
		normcheck+=COAL_DELR*4.0*PI*r*r*psi_coal[ir]*psi_coal[ir];
	}
	if(fabs(normcheck-1.0)>1.0E-5){
		printf("normcheck=%g\n",normcheck);
		exit(1);
	}
}

void CF::WriteCFs(string filename){
	FILE *fptr=fopen(filename.c_str(),"w");
	int iq;
	fprintf(fptr,"----- CF ------, norm_qinv[0]=%g\n",norm_qinv[0]);
	fprintf(fptr,"q(MeV/c) CF(qinv) CF(qout) CF(side) CF(qlong)\n");
	for(iq=0;iq<NQ;iq++){
		fprintf(fptr,"%7.3f %8.5f %8.5f %8.5f %8.5f ",(iq+0.5)*DELQ,cf_qinv[iq],cf_qout[iq],cf_qside[iq],cf_qlong[iq]);
		fprintf(fptr,"%15.8e %15.8e %15.8e %15.8e\n",norm_qinv[iq],norm_qout[iq],norm_qside[iq],norm_qlong[iq]);
	}
	fclose(fptr);
}

void CF::ReadCFs(string filename){
	char dummy[120];
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,120,fptr);
	fgets(dummy,120,fptr);
	int iq;
	double qinv;
	for(iq=0;iq<NQ;iq++){
		fscanf(fptr,"%lf %lf %lf %lf %lf",&qinv,&cf_qinv[iq],&cf_qout[iq],&cf_qside[iq],&cf_qlong[iq]);
		fscanf(fptr,"%lf %lf %lf %lf",&norm_qinv[iq],&norm_qout[iq],&norm_qside[iq],&norm_qlong[iq]);
	}
	fclose(fptr);
}

#endif