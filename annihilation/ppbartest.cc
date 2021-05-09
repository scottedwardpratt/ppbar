#include <coral.h>
#include "hbt_bes.h"

using namespace std;

int main(int argc,char *argv[]){
	string parsfilename;
	if(argc!=2){
		printf("usage: ppbartest parsfilename\n");
		exit(1);
	}
	parsfilename=argv[1];
	double q,r,ctheta,x,y,z,Rx,Ry,Rz,psi2,delq;
	vector<double>psi2bar;
	unsigned long long int imc,nmc;
	int iq,nqmax;
	CWaveFunction *wf;
	//wf=new CWaveFunction_ppbar_nocoulomb(parsfilename);
	wf=new CWaveFunction_optical(parsfilename);
	printf("wf created\n");
	nmc=wf->parameters.getI("NMC",1000);
	CRandy *randy=new CRandy(time(NULL));
	delq=wf->GetDELQ();
	nqmax=wf->GetNQMAX();
	psi2bar.resize(nqmax);
	for(iq=0;iq<nqmax;iq++){
		psi2bar[iq]=0.0;
	}
	Rx=Ry=Rz=6.0;
	for(imc=0;imc<nmc;imc++){
		x=Rx*randy->ran_gauss();
		y=Ry*randy->ran_gauss();
		z=Rz*randy->ran_gauss();
		r=sqrt(x*x+y*y+z*z);
		ctheta=1.0-2.0*randy->ran();
		for(iq=0;iq<nqmax;iq++){
			psi2=wf->GetPsiSquared(iq,r,ctheta);
			psi2bar[iq]+=psi2;
		}
	}
	char filename[100];
	sprintf(filename,"CF_Rx%gRy%gRz%g.txt",Rx,Ry,Rz);
	FILE *fptr=fopen(filename,"w");
	for(iq=0;iq<nqmax;iq++){
		q=(iq+0.5)*delq;
		printf("%5.1f %g\n",q,psi2bar[iq]/double(nmc));
		fprintf(fptr,"%5.1f %g\n",q,psi2bar[iq]/double(nmc));
	}
	fclose(fptr);
	delete wf;
	return 0;
}
