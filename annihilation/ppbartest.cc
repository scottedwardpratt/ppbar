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
	char filename[100];
	FILE *fptr;
	double q,r,ctheta,x,y,z,Rx,Ry,Rz,psi2,delq;
	vector<double>psi2bar;
	unsigned long long int imc,nmc;
	int iq,nqmax;
	CWaveFunction_optical *wf;
	wf=new CWaveFunction_optical(parsfilename);
	nmc=wf->parameters.getI("NMC",1000);
	CRandy *randy=new CRandy(time(NULL));
	delq=wf->GetDELQ();
	nqmax=wf->GetNQMAX();
	psi2bar.resize(nqmax);
	Rx=Ry=Rz=6.0;
	sprintf(filename,"CF_Rx%gRy%gRz%g.txt",Rx,Ry,Rz);
	fptr=fopen(filename,"w");
	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		psi2bar[iq]=0.0;
		for(imc=0;imc<nmc;imc++){
			x=Rx*randy->ran_gauss();
			y=Ry*randy->ran_gauss();
			z=Rz*randy->ran_gauss();
			r=sqrt(x*x+y*y+z*z);
			ctheta=1.0-2.0*randy->ran();
			psi2=wf->GetPsiSquared(q,r,ctheta);
			psi2bar[iq]+=psi2;
		}
		printf("%5.1f %8.5f %7.2f mb\n",q,psi2bar[iq]/double(nmc),10.0*wf->sigma_annihilation[iq]);
		fprintf(fptr,"%5.1f %8.5f %7.2f mb\n",q,psi2bar[iq]/double(nmc),10.0*wf->sigma_annihilation[iq]);
	}
	fclose(fptr);
	delete wf;
	return 0;
}
