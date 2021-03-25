#include "coral.h"
#include "commondefs.h"
#include "hbt_bes.h"
#include "wavefunction.h"
#include "constants.h"
#include "sf.h"

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

int main(){
	CRandy *randy=new CRandy(-1234);
	double root2=sqrt(2.0);
	double x[3];
	x[0]=root2*3*randy->ran_gauss();
	x[1]=root2*4*randy->ran_gauss();
	x[2]=root2*5*randy->ran_gauss();
	int iq;
	double phi,ctheta,stheta,qx,qy,qz,q;
	double r,psisquared,psi,ctheta_qr;
	string filename="parameters/parameters_pp.txt";
	CWaveFunction *wf=new CWaveFunction_generic(filename, -1, ProtonMass, ProtonMass, 0.5);
	r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	
	
	for(iq=0; iq<CF::NSAMPLE_THETAPHI; iq++){
		phi=2.0*PI*randy->ran();
		ctheta=-1.0+2.0*randy->ran();
		stheta=sqrt(1.0-ctheta*ctheta);
		qz=ctheta;
		qx=stheta*cos(phi);
		qy=stheta*sin(phi);
		q=sqrt(qx*qx+qy*qy+qz*qz);
		ctheta_qr=(qx*x[0]+qy*x[1]+qz*x[2])/r;
		
		psisquared=wf->GetPsiSquared(q,r,ctheta_qr);
		psi=wf->GetPsiSquared(q,r,-ctheta_qr);
		if(psisquared != psi){
			printf("Psi error");
			exit (1);
		}
	}
	
	return 0;
}

