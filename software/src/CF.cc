#ifndef __HBT_BES_CALC_CF_CC__
#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

int CF::NQ=50;
double CF::DELQ=2.0;
CHBT_BES *CF::hbt=NULL;
double CF::OUTSIDELONG_DIRECTION_CUT=0.9;
double CF::OUTSIDELONG_Q_CUT=10.0;
bool CF::USE_OUTSIDELONG_Q_CUT=true;
bool CF::USE_OUTSIDELONG_DIRECTION_CUT=false;
CRandy* CF::randy=NULL;

CF::CF(){
	cf_qinv.resize(NQ);
	cf_qout.resize(NQ);
	cf_qside.resize(NQ);
	cf_qlong.resize(NQ);
	ThetaPhiDist.resize(10);
	for(int ictheta=0;ictheta<10;ictheta++)
		ThetaPhiDist[ictheta].resize(18);
	Reset();
}

void CF::Reset(){
	nsample_qinv=nsample_qout=nsample_qside=nsample_qlong=0;
	for(int iq=0;iq<NQ;iq++){
		cf_qinv[iq]=cf_qout[iq]=cf_qside[iq]=cf_qlong[iq]=0.0;
	}
	for(int ictheta=0;ictheta<10;ictheta++)
		for(int iphi=0;iphi<18;iphi++)
			ThetaPhiDist[ictheta][iphi]=0.0;
}

void CF::Increment(CHBT_Part *parta,CHBT_Part *partb){
	int iq,ithetaphi,iphi,ictheta;
	const int Nthetaphi=10;
	double phi,ctheta,stheta,qx,qy,qz,qinv;
	double r,psisquared,ctheta_qr;
	vector<double> x(4,0.0);
	CalcXR(parta,partb,x,r);
	//x[1]=7.0*randy->ran_gauss();
	//x[2]=4.0*randy->ran_gauss();
	//x[3]=6.0*randy->ran_gauss();
	//r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
	if(r==r && r!=0.0){
		// Increment correlation function
		for(ithetaphi=0;ithetaphi<Nthetaphi;ithetaphi++){
			phi=2.0*PI*randy->ran();
			ctheta=-1.0+2.0*randy->ran();
			stheta=sqrt(1.0-ctheta*ctheta);
			qz=ctheta;
			qx=stheta*cos(phi);
			qy=stheta*sin(phi);
			qinv=(iq+0.5)*DELQ;
			ctheta_qr=(qx*x[1]+qy*x[2]+qz*x[3])/r;		
			for(iq=0;iq<NQ;iq++){
				psisquared=wf->CalcPsiSquared(iq,r,ctheta_qr);
				//psisquared=1.0;
				if(psisquared!=psisquared){
					parta->Print();
					partb->Print();
					printf("---------- psisquared = %g ------------\n",psisquared);
					exit(1);
				}
				cf_qinv[iq]+=psisquared;
				nsample_qinv+=1;
				if(USE_OUTSIDELONG_DIRECTION_CUT){
					if(fabs(qx)>OUTSIDELONG_DIRECTION_CUT){
						cf_qout[iq]+=psisquared;
						nsample_qout+=1;
					}
					if(fabs(qy)>OUTSIDELONG_DIRECTION_CUT){
						cf_qside[iq]+=psisquared;
						nsample_qside+=1;
					}
					if(fabs(qz)>OUTSIDELONG_DIRECTION_CUT){
						cf_qlong[iq]+=psisquared;
						nsample_qlong+=1;
					}
				}
				if(USE_OUTSIDELONG_Q_CUT){
					double qperp;
					qperp=qinv*sqrt(1.0-qx*qx);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qout[iq]+=psisquared;
						nsample_qout+=1;
					}
					qperp=qinv*sqrt(1.0-qy*qy);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qside[iq]+=psisquared;
						nsample_qside+=1;
					}
					qperp=qinv*sqrt(1.0-qz*qz);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qlong[iq]+=psisquared;
						nsample_qlong+=1;
					}
				}
			}
		}
		// Increment ThetaPhiDist
		if(r<20.0){
			ctheta=x[3]/r;
			phi=atan2(x[2],x[1]);
			if(phi<0)
				phi+=2.0*PI;
			if(phi>PI){
				phi=phi-PI;
				ctheta=-ctheta;
			}
			ictheta=lrint(floor(5.0*(1.0+ctheta)));
			iphi=lrint(floor(phi*18.0/PI));
			if(ictheta<0 || ictheta>=int(ThetaPhiDist.size())){
				printf("ictheta too big=%d\n",ictheta);
				exit(1);
			}
			if(iphi<0 || iphi>=int(ThetaPhiDist[0].size())){
				printf("iphi too big=%d\n",iphi);
				exit(1);
			}
			ThetaPhiDist[ictheta][iphi]+=1.0;
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
	if(nsample_qinv>0){
		for(iq=0;iq<NQ;iq++){
			cf_qinv[iq]=NQ*cf_qinv[iq]/double(nsample_qinv);
		}
	}
	if(nsample_qout>0){
		for(iq=0;iq<NQ;iq++){
			cf_qout[iq]=NQ*cf_qout[iq]/double(nsample_qout);
		}
	}
	if(nsample_qside>0){
		for(iq=0;iq<NQ;iq++){
			cf_qside[iq]=NQ*cf_qside[iq]/double(nsample_qside);
		}
	}
	if(nsample_qlong>0){
		for(iq=0;iq<NQ;iq++){
			cf_qlong[iq]=NQ*cf_qlong[iq]/double(nsample_qlong);
		}
	}
}

void CF::Print(){
	int iq;
	printf("----- CF ------, nsample_qinv=%d\n",nsample_qinv);
	printf("q(MeV/c) CF(qinv) CF(qout) CF(side) CF(qlong)\n");
	for(iq=0;iq<NQ;iq++){
		printf("%7.3f %8.5f %8.5f %8.5f %8.5f\n",(iq+0.5)*DELQ,cf_qinv[iq],cf_qout[iq],cf_qside[iq],cf_qlong[iq]);
	}
}

void CF::WriteCFs(string filename){
	FILE *fptr=fopen(filename.c_str(),"w");
	int iq;
	fprintf(fptr,"----- CF ------, nsample_qinv=%d\n",nsample_qinv);
	fprintf(fptr,"q(MeV/c) CF(qinv) CF(qout) CF(side) CF(qlong)\n");
	for(iq=0;iq<NQ;iq++){
		fprintf(fptr,"%7.3f %8.5f %8.5f %8.5f %8.5f\n",(iq+0.5)*DELQ,cf_qinv[iq],cf_qout[iq],cf_qside[iq],cf_qlong[iq]);
	}
	fclose(fptr);
}

void CF::WriteThetaPhiDist(string filename){
	FILE *fptr=fopen(filename.c_str(),"w");
	int ictheta,iphi;
	double norm=0.0;
	for(ictheta=0;ictheta<10;ictheta++){
		for(iphi=0;iphi<18;iphi++){
			norm+=ThetaPhiDist[ictheta][iphi];
		}
	}
	norm=norm/180.0;
	
	fprintf(fptr,"cos(theta)|          ------  PHI (degrees) ------\n");
	for(iphi=0;iphi<18;iphi++){
		fprintf(fptr,"  %5.1f   ",(iphi+0.5)*10);
	}
	fprintf(fptr,"\n");
	for(ictheta=0;ictheta<10;ictheta++){
		fprintf(fptr,"%5.1f",-1.0+2.0*(ictheta+0.5)/10.0);
		for(iphi=0;iphi<18;iphi++){
			fprintf(fptr," %8.5f ",ThetaPhiDist[ictheta][iphi]/norm);
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
}


#endif