#ifndef __HBT_SPECTRA_CC__
#define __HBT_SPECTRA_CC__
#include <cstring>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

void CHBT_BES::CalcWriteSpectra(){
	double COAL_SPIN_FACTOR=parmap->getD("COAL_SPIN_FACTOR",0.75);
	double Pt,uperp,D3PoverE,D3PoverEa,D3PoverEb,E,phi,delN,uperp1,uperp0,mu=MASSA*MASSB/(MASSA+MASSB);;
	double Ebar=0.0,Ptbar=0.0,EtPtnorm=0.0,V2=0.0;
	int iuperp,irap,iphi;
	CF *cfptr;
	double cfactor;
	//D3r=4.0*PI*pow(CF::Rcoalescence,3)/3.0; // Coalescence volume in coordinate space
	vector<double> coalspectra,aspectra,bspectra,Rhbt,B2;
	coalspectra.resize(NUPERP);
	aspectra.resize(NUPERP);
	bspectra.resize(NUPERP);
	Rhbt.resize(NUPERP);
	B2.resize(NUPERP);
	for(iuperp=0;iuperp<NUPERP;iuperp++){
		coalspectra[iuperp]=aspectra[iuperp]=bspectra[iuperp]=0.0;
	}
	
	long long int ntot=0;
	double ncoalraw=0,ncoal=0.0;
	cfactor=COAL_SPIN_FACTOR*pow(2.0*PI*HBARC,3)/(double(NEVENTS)*double(NEVENTS));
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			phi=DELPHI*(iphi+0.5);
			ntot+=partmap[irap][iphi].size();
			for(iuperp=0;iuperp<NUPERP;iuperp++){
				cfptr=CFArray[irap][iphi][iuperp];
				uperp=(iuperp+0.5)*DELUPERP;
				Pt=(MASSA+MASSB)*uperp;
				E=sqrt(Pt*Pt+(MASSA+MASSB)*(MASSA+MASSB));
				delN=cfactor*cfptr->ncoalescence;
				ncoal+=delN;
				ncoalraw+=cfptr->ncoalescence;
				if(IDA==IDB)
					coalspectra[iuperp]+=2.0*delN;
				else
					coalspectra[iuperp]+=delN;
				EtPtnorm+=delN;
				Ebar+=E*delN;
				Ptbar+=Pt*delN;
				V2+=cos(2.0*phi)*delN;
				aspectra[iuperp]+=Na[irap][iphi][iuperp];
				bspectra[iuperp]+=Nb[irap][iphi][iuperp];
			}
		}
	}
	printf("ncoal=%g, ncoalraw=%g\n",ncoal,ncoalraw);
	for(iuperp=0;iuperp<NUPERP;iuperp++){
		Pt=(iuperp+0.5)*DELPT;
		D3PoverE=2.0*PI*NRAP*DELRAP*Pt*DELPT;
		coalspectra[iuperp]=coalspectra[iuperp]/D3PoverE;
		printf("%6.2f %g\n",Pt,coalspectra[iuperp]*1.0E6);
	}
	printf("For coalescence, <Et>=%g, <Pt>=%g, V2=%g\n",Ebar/EtPtnorm,Ptbar/EtPtnorm,V2/EtPtnorm);
	printf("ntot=%lld\n",ntot);
	mu=MASSA*MASSB/(MASSA+MASSB);
	cfactor=pow(PI*HBARC*HBARC,1.5)*COAL_SPIN_FACTOR/mu;
	char command[150],filename[150];
	sprintf(command,"mkdir -p %s/spectra",RESULTS_DIR.c_str());
	system(command);
	sprintf(filename,"%s/spectra/spectra_%d.txt",RESULTS_DIR.c_str(),IDA);
	FILE *fptra=fopen(filename,"w");
	sprintf(filename,"%s/spectra/spectra_%d.txt",RESULTS_DIR.c_str(),IDB);
	FILE *fptrb=fopen(filename,"w");
	sprintf(filename,"%s/spectra/spectra_coalescence.txt",RESULTS_DIR.c_str());
	FILE *fptrc=fopen(filename,"w");
	fprintf(fptrc," Pt     Rhbt        B2        spectra\n");
	for(iuperp=0;iuperp<NUPERP;iuperp++){
		uperp0=iuperp*DELUPERP;
		uperp1=(iuperp+1)*DELUPERP;
		D3PoverEa=DELRAP*NRAP*MASSA*MASSA*PI*(uperp1*uperp1-uperp0*uperp0);
		D3PoverEb=DELRAP*NRAP*MASSB*MASSB*PI*(uperp1*uperp1-uperp0*uperp0);
		aspectra[iuperp]=aspectra[iuperp]/(NEVENTS*D3PoverEa);
		bspectra[iuperp]=bspectra[iuperp]/(NEVENTS*D3PoverEb);
		Rhbt[iuperp]=cfactor*aspectra[iuperp]*bspectra[iuperp]/coalspectra[iuperp];
		B2[iuperp]=coalspectra[iuperp]/(aspectra[iuperp]*bspectra[iuperp]);
		Rhbt[iuperp]=pow(Rhbt[iuperp],1.0/3.0);
		Pt=MASSA*(iuperp+0.5)*DELUPERP;
		fprintf(fptra,"%7.1f %15.7e\n",Pt,aspectra[iuperp]*1.0E6);
		Pt=MASSB*(iuperp+0.5)*DELUPERP;
		fprintf(fptrb,"%7.1f %15.7e\n",Pt,bspectra[iuperp]*1.0E6);
		Pt=(MASSA+MASSB)*(iuperp+0.5)*DELUPERP;
		printf("%7.1f %7.4f %15.7e %15.7e\n",Pt,Rhbt[iuperp],B2[iuperp]/1.0E6,coalspectra[iuperp]*1.0E6);
		fprintf(fptrc,"%7.1f %7.4f %15.7e %15.7e\n",Pt,Rhbt[iuperp],B2[iuperp]/1.0E6,coalspectra[iuperp]*1.0E6);
	}
	fclose(fptra);
	fclose(fptrb);
	fclose(fptrc);
}

#endif
