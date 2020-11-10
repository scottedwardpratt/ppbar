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
double CF::Dxyz=0.5;
int CF::Nxyz=100;
CRandy* CF::randy=NULL;

CF::CF(){
	cf_qinv.resize(NQ);
	cf_qout.resize(NQ);
	cf_qside.resize(NQ);
	cf_qlong.resize(NQ);
	norm_qinv.resize(NQ);
	norm_qout.resize(NQ);
	norm_qside.resize(NQ);
	norm_qlong.resize(NQ);
	ThetaPhiDist.resize(10);
	source_out.resize(Nxyz);
	source_side.resize(Nxyz);
	source_long.resize(Nxyz);
	for(int ictheta=0;ictheta<10;ictheta++)
		ThetaPhiDist[ictheta].resize(18);
	Reset();
}

void CF::Reset(){
	nincrement=0;
	for(int iq=0;iq<NQ;iq++){
		cf_qinv[iq]=cf_qout[iq]=cf_qside[iq]=cf_qlong[iq]=0.0;
		norm_qinv[iq]=norm_qout[iq]=norm_qside[iq]=norm_qlong[iq]=0;
	}
	for(int ictheta=0;ictheta<10;ictheta++)
		for(int iphi=0;iphi<18;iphi++)
			ThetaPhiDist[ictheta][iphi]=0.0;
	for(int ixyz=0;ixyz<Nxyz;ixyz++){
		source_out[ixyz]=source_side[ixyz]=source_long[ixyz]=0.0;
	}
}

void CF::Increment(CHBT_Part *parta,CHBT_Part *partb){
	int iq,ithetaphi,iphi,ictheta,ixyz;
	const int Nthetaphi=10;
	double phi,ctheta,stheta,qx,qy,qz,qinv;
	double r,psisquared,ctheta_qr;
	vector<double> x(4,0.0);
	CalcXR(parta,partb,x,r);
	//x[1]=7.0*randy->ran_gauss();
	//x[2]=5.0*randy->ran_gauss();
	//x[3]=5.0*randy->ran_gauss();
	//r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
	nincrement+=1;
	ixyz=fabs(x[1])/Dxyz;
	if(ixyz<Nxyz)
		source_out[ixyz]+=1.0;
	ixyz=fabs(x[2])/Dxyz;
	if(ixyz<Nxyz)
		source_side[ixyz]+=1.0;
	ixyz=fabs(x[3])/Dxyz;
	if(ixyz<Nxyz)
		source_long[ixyz]+=1.0;
	if(r==r && r!=0.0){
		// Increment correlation function
		for(ithetaphi=0;ithetaphi<Nthetaphi;ithetaphi++){
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
				cf_qinv[iq]+=psisquared;
				norm_qinv[iq]+=1;
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
						cf_qout[iq]+=psisquared;
						norm_qout[iq]+=1;
					}
					qperp=qinv*sqrt(1.0-qy*qy);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qside[iq]+=psisquared;
						norm_qside[iq]+=1;
					}
					qperp=qinv*sqrt(1.0-qz*qz);
					if(qperp<OUTSIDELONG_Q_CUT){
						cf_qlong[iq]+=psisquared;
						norm_qlong[iq]+=1;
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

void CF::Increment(vector<double> &x){

	int iq,ithetaphi,iphi,ictheta;
	const int Nthetaphi=10;
	double phi,ctheta,stheta,qx,qy,qz,qinv,qperp;
	double r,psisquared,ctheta_qr;
	r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
	nincrement+=1;

	for(ithetaphi=0;ithetaphi<Nthetaphi;ithetaphi++){
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
			cf_qinv[iq]+=psisquared;
			norm_qinv[iq]+=1;
			if(USE_OUTSIDELONG_Q_CUT){
				qperp=qinv*sqrt(1.0-qx*qx);
				if(qperp<OUTSIDELONG_Q_CUT){
					cf_qout[iq]+=psisquared;
					norm_qout[iq]+=1;
				}
				qperp=qinv*sqrt(1.0-qy*qy);
				if(qperp<OUTSIDELONG_Q_CUT){
					cf_qside[iq]+=psisquared;
					norm_qside[iq]+=1;
				}
				qperp=qinv*sqrt(1.0-qz*qz);
				if(qperp<OUTSIDELONG_Q_CUT){
					cf_qlong[iq]+=psisquared;
					norm_qlong[iq]+=1;
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
			cf_qinv[iq]=cf_qinv[iq]/double(norm_qinv[iq]);
		}
		else
			cf_qinv[iq]=0;
		if(norm_qout[iq]>0){
			cf_qout[iq]=cf_qout[iq]/double(norm_qout[iq]);
		}
		else
			cf_qout[iq]=0;
		if(norm_qside[iq]>0){
			cf_qside[iq]=cf_qside[iq]/double(norm_qside[iq]);
		}
		else
			cf_qside[iq]=0;
		if(norm_qlong[iq]>0){
			cf_qlong[iq]=cf_qlong[iq]/double(norm_qlong[iq]);
		}
		else
			cf_qlong[iq]=0;
	}
}

void CF::Print(){
	int iq,ixyz;
	printf("#---- CF ------, norm_qinv=%lld\n",norm_qinv[0]);
	printf("q(MeV/c) CF(qinv) CF(qout) CF(side) CF(qlong)\n");
	for(iq=0;iq<NQ;iq++){
		printf("%7.3f %8.5f %8.5f %8.5f %8.5f\n",(iq+0.5)*DELQ,cf_qinv[iq],cf_qout[iq],cf_qside[iq],cf_qlong[iq]);
	}
	printf("nincrement=%lld\n",nincrement);
	printf("# xyz     source_out     source_side   source_long\n");
	for(ixyz=0;ixyz<Nxyz;ixyz++){
		printf("%5.2f %10.4e %10.4e %10.4e\n",(ixyz+0.5)*Dxyz,source_out[ixyz],source_side[ixyz],source_long[ixyz]);
	}
}

void CF::WriteCFs(string filename){
	FILE *fptr=fopen(filename.c_str(),"w");
	int iq;
	fprintf(fptr,"----- CF ------, norm_qinv=%lld\n",norm_qinv[0]);
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