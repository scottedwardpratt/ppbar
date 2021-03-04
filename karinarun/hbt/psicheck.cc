//#ifndef __HBT_BES_CALC_CF_CC__
#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

int CF::NQ=40;
double CF::DELQ=2.0;
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
int CF::NSAMPLE_THETAPHI=4;

int main(int argc, char *argv[]){
	double root2=sqrt(2.0);
	double x[3]= {root2*3*randy->ran_gauss(),root2*4*randy->ran_gauss(),root2*5*randy->ran_gauss()};
	int iq,ithetaphi;
	double phi,ctheta,stheta,qx,qy,qz,qinv,qperp, cosgamma;
	double r,psisquared,psi,ctheta_qr;
	r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
	phi=2.0*PI*randy->ran();
	ctheta=-1.0+2.0*randy->ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	qz=ctheta;
	qx=stheta*cos(phi);
	qy=stheta*sin(phi);
	ctheta_qr=(qx*x[1]+qy*x[2]+qz*x[3])/r;
	
	for(int iq=0; iq<1000; iq++){
		psisquared=wf->CalcPsiSquared(iq,r,ctheta_qr);
		psi=wf->CalcPsiSquared(-iq,-r,ctheta_qr);
		if(psisquared != psi){
			printf("Psi error");
			exit;
		}
	}
	
	return 0;
}
