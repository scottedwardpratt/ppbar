#ifndef __HBT_BES_CC__
#define __HBT_BES_CC__
#include <cstring>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

CHBT_BES::CHBT_BES(string parsfilename){
	int irap,iphi,iuperp;
	double D3q,mu,uperp;
	GITHOME_MSU=getenv("GITHOME_MSU");
	printf("GITHOME_MSU=%s\n",GITHOME_MSU.c_str());
	parmap=new CparameterMap();
	parmap->ReadParsFromFile(parsfilename);
	TAU_COMPARE=parmap->getD("TAU_COMPARE",12.0);
	GAUSS=parmap->getB("GAUSS",false);
	CF::NQ=parmap->getI("NQMAX",100);
	CF::DELQ=parmap->getD("DELQ",2.0);
	CF::OUTSIDELONG_DIRECTION_CUT=parmap->getD("OUTSIDELONG_DIRECTION_CUT",0.9);
	CF::OUTSIDELONG_Q_CUT=parmap->getD("OUTSIDELONG_Q_CUT",10.0);
	CF::USE_OUTSIDELONG_DIRECTION_CUT=parmap->getB("USE_OUTSIDELONG_DIRECTION_CUT",false);
	CF::USE_OUTSIDELONG_Q_CUT=parmap->getB("USE_OUTSIDELONG_Q_CUT",true);
	CF::Nxyz=100;
	CF::Dxyz=0.5;
	CF::Rcoalescence=parmap->getD("RCOALESCENCE",2.0);
	CF::NSAMPLE_THETAPHI=parmap->getI("NSAMPLE_THETAPHI",4);
	CF::COAL_USE_WF=parmap->getB("COAL_USE_WF",false);
	RESULTS_DIR=parmap->getS("RESULTS_DIR","results");
	NPHI=parmap->getI("NPHI",16);
	DELPHI=2.0*PI/double(NPHI);
	NRAP=parmap->getI("NRAP",3);
	DELRAP=parmap->getD("DELRAP",0.2);
	NUPERP=parmap->getI("NUPERP",150);
	DELPT=parmap->getD("DELPT",50.0);
	IDA=parmap->getI("IDA",211);
	IDB=parmap->getI("IDB",211);
	MASSA=parmap->getD("MASSA",-12345.0);
	MASSB=parmap->getD("MASSB",-12345.0);
	QAB=parmap->getD("QAB",1);
	NEVENTS_MAX=parmap->getI("NEVENTS_MAX",10);
	NMC=parmap->getI("NMC",100000);
	NTRY=0;
	TAU_COMPARE=25.0;
	INPUT_OSCAR_BASE_DIRECTORY=GITHOME_MSU+"/"+parmap->getS("INPUT_OSCAR_BASE_DIRECTORY","crap");
	INPUT_OSCAR_NRUNS=parmap->getI("OSCAR_NRUNS",24);
	//RANSEED=parmap->getD("RANSEED",-12345);
	RANSEED=-time(NULL);
	UPERPTEST=parmap->getD("UPERPTEST",0.02);
	randy=new CRandy(RANSEED);
	CF::randy=randy;
	CF::hbt=this;
	Na.resize(NRAP);
	Nb.resize(NRAP);
	for(irap=0;irap<NRAP;irap++){
		Na[irap].resize(NPHI);
		Nb[irap].resize(NPHI);
		for(iphi=0;iphi<NPHI;iphi++){
			Na[irap][iphi].resize(NUPERP);
			Nb[irap][iphi].resize(NUPERP);
			for(iuperp=0;iuperp<NUPERP;iuperp++){
				Na[irap][iphi][iuperp]=0;
				Nb[irap][iphi][iuperp]=0;
			}
		}
	}
	
	if((IDA==2212 && IDB==2212) || (IDA==-2212 && IDB==-2212)){
		wf=new CWaveFunction_pp_schrod(parsfilename);
		MASSA=MASSB=938.272;
	}
	else if((IDA==211 && IDB==211) || (IDA==-211 && IDB==-211)){
		wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);
		MASSA=MASSB=139.570;
	}
	else if((IDA==211 && IDB==-211) || (IDA==-211 && IDB==211)){
		wf=new CWaveFunction_pipluspiminus_sqwell(parsfilename);
		MASSA=MASSB=139.570;
	}
	else if((IDA==321 && IDB==321) || (IDA==-321 && IDB==-321)){
		MASSA=MASSB=493.677;
		wf=new CWaveFunction_generic(parsfilename,QAB,MASSA,MASSB,1.0);
	}
	else if((IDA==321 && IDB==-321) || (IDA==-321 && IDB==321)){
		MASSA=MASSB=493.677;
		wf=new CWaveFunction_generic(parsfilename,QAB,MASSA,MASSB,0.5);
	}
	else if((IDA==2212 && IDB==-2212) || (IDA==-2212 && IDB==2212)){
		MASSA=MASSB=938.272;
		wf=new CWaveFunction_generic(parsfilename,QAB,MASSA,MASSB,0.5);
	}
	else{
		printf("fatal: IDs not recognized, IDA=%d, IDB=%d\n",IDA,IDB);
		exit(1);
	}
	partmap.resize(NRAP);
	for(irap=0;irap<NRAP;irap++){
		partmap[irap].resize(NPHI);
	}
	DELUPERP=DELPT/(MASSA+MASSB);
	
	if(!GAUSS){
		CFArray.resize(NRAP);
		for(irap=0;irap<NRAP;irap++){
			CFArray[irap].resize(NPHI);
			for(iphi=0;iphi<NPHI;iphi++){
				CFArray[irap][iphi].resize(NUPERP);
				for(iuperp=0;iuperp<NUPERP;iuperp++){
					CFArray[irap][iphi][iuperp]=new CF();
					CFArray[irap][iphi][iuperp]->Reset();
					CFArray[irap][iphi][iuperp]->wf=wf;
					uperp=(iuperp+0.5)*DELUPERP;
					mu=MASSA*MASSB/(MASSA+MASSB);
					D3q=2.0*mu*DELRAP;
					D3q*=2.0*mu*uperp*DELPHI;
					D3q*=2.0*mu*UPERPTEST;
					CFArray[irap][iphi][iuperp]->D3q=D3q;
				}
			}
		}
		cfbar=new CF();
		cfbar->Reset();
		cfbar->wf=wf;
		if(CF::COAL_USE_WF)
			cfbar->CalcCoalWF();
	}
	else{
		cfgauss=new CF();
		cfgauss->nincrement=0;
		cfgauss->Reset();
		cfgauss->wf=wf;
		cfgauss->D3q=1.0;
	}
}

void CHBT_BES::AddPart(int &IDread,vector<double> &p,vector<double> &x,double mass){
	CHBT_Part *newpart;
	double rap0,phi0,uperp,pt,y,z,cphi,sphi,gg,ggv;
	int irap,iphi,iuperp;
	rap0=atanh(p[3]/p[0]);
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	uperp=pt/mass;	
	phi0=atan2(p[2],p[1]);
	GetIrapIphiIuperp(rap0,phi0,uperp,irap,iphi,iuperp);
	for(int alpha=0;alpha<4;alpha++){
		p[alpha]*=1000.0;
	}
	mass*=1000.0;
	
	if(irap>=0 && irap<NRAP && iuperp<NUPERP){
		newpart=new CHBT_Part();
		newpart->ID=IDread;
		cphi=cos(phi0);
		sphi=sin(phi0);
		y=x[2];
		x[2]=x[2]*cphi-x[1]*sphi;
		x[1]=x[1]*cphi+y*sphi;
		p[2]=0.0;
		p[1]=pt;
		gg=cosh(rap0);
		ggv=sinh(rap0);
		p[0]=p[0]*gg-ggv*p[3];
		p[3]=0.0;
		z=x[3];
		x[3]=gg*x[3]-ggv*x[0];
		x[0]=gg*x[0]-ggv*z;
		x[1]=x[1]-(p[1]/p[0])*(x[0]-TAU_COMPARE);
		x[0]=TAU_COMPARE;
		for(int alpha=0;alpha<4;alpha++){
			newpart->p[alpha]=p[alpha];
			newpart->x[alpha]=x[alpha];
		}
		newpart->mass=sqrt(p[0]*p[0]-pt*pt-p[3]*p[3]);
		newpart->pt=pt;
		newpart->uperp=uperp;
		newpart->rap0=rap0;
		newpart->phi0=phi0;
		partmap[irap][iphi].insert(pair<double,CHBT_Part*>(uperp,newpart));
		if(newpart->ID==IDA)
			Na[irap][iphi][iuperp]+=1;
		if(newpart->ID==IDB)
			Nb[irap][iphi][iuperp]+=1;
	}
}

double CHBT_BES::Getqinv(vector<double> &pa,vector<double> &pb){
	vector<double> P,q;
	P.resize(4); q.resize(4);
	double M2,Pdotq,QINV2;
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=pa[alpha]+pb[alpha];
		q[alpha]=pa[alpha]-pb[alpha];
	}
	M2=P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3];
	Pdotq=P[0]*q[0]-P[1]*q[1]-P[2]*q[2]-P[3]*q[3];
	for(alpha=0;alpha<4;alpha++)
		q[alpha]-=Pdotq*P[alpha]/M2;
	QINV2=q[1]*q[1]+q[2]*q[2]+q[3]*q[3]-q[0]*q[0];
	return 0.5*sqrt(QINV2);
}

CF* CHBT_BES::GetCF(CHBT_Part *partaa,CHBT_Part *partbb){
	int irap,iphi,iuperp;
	double phibar,delphi;
	delphi=partaa->phi0-partbb->phi0;
	phibar=0.5*(partaa->phi0+partbb->phi0);
	if(fabs(delphi)>PI){
		if(phibar>0.0)
			phibar=phibar-PI;
		if(phibar<0.0)				
			phibar=phibar+PI;
	}
	GetIrapIphiIuperp(0.5*(partaa->rap0+partbb->rap0),phibar,0.5*(partaa->uperp+partbb->uperp),
	irap,iphi,iuperp);
	if(irap<NRAP && iuperp<NUPERP){
		return CFArray[irap][iphi][iuperp];
	}
	else
		return NULL;
}

void CHBT_BES::GetIrapIphiIuperp(double rap,double phi,double uperp,int &irap,int &iphi,int &iuperp){
	irap=floorl((0.5*NRAP*DELRAP+rap)/DELRAP);
	iphi=floorl((PI+phi)/DELPHI);
	iuperp=floorl(uperp/DELUPERP);
}

void CHBT_BES::AverageCF(){
	CF *cfptr;
	int iq;
	cfbar->Reset();
	for(int irap=0;irap<NRAP;irap++){
		for(int iphi=0;iphi<NPHI;iphi++){
			for(int iuperp=0;iuperp<NUPERP;iuperp++){
				cfptr=CFArray[irap][iphi][iuperp];
				cfbar->nincrement+=cfptr->nincrement;
				if(cfptr->norm_qinv[0]>1.0E-20){
					for(iq=0;iq<CF::NQ;iq++){
						cfbar->norm_qinv[iq]+=cfptr->norm_qinv[iq];
						cfbar->norm_qout[iq]+=cfptr->norm_qout[iq];
						cfbar->norm_qside[iq]+=cfptr->norm_qside[iq];
						cfbar->norm_qlong[iq]+=cfptr->norm_qlong[iq];
						cfbar->cf_qinv[iq]+=cfptr->cf_qinv[iq]*cfptr->norm_qinv[iq];
						cfbar->cf_qout[iq]+=cfptr->cf_qout[iq]*cfptr->norm_qout[iq];
						cfbar->cf_qside[iq]+=cfptr->cf_qside[iq]*cfptr->norm_qside[iq];
						cfbar->cf_qlong[iq]+=cfptr->cf_qlong[iq]*cfptr->norm_qlong[iq];
						
						if(cfptr->cf_qinv[iq] != cfptr->cf_qinv[iq]){
							printf("irap=%d, iphi=%d, iuperp=%d\n",irap,iphi,iuperp);
							CFArray[irap][iphi][iuperp]->Print();
							exit(1);
						}
					}
					for(int ictheta=0;ictheta<10;ictheta++){
						for(int iphir=0;iphir<18;iphir++){
							cfbar->ThetaPhiDist[ictheta][iphir]+=cfptr->ThetaPhiDist[ictheta][iphir]*cfptr->nincrement;
						}
					}
				}
				for(int ixyz=0;ixyz<CF::Nxyz;ixyz++){
					cfbar->source_out[ixyz]+=cfptr->source_out[ixyz];
					cfbar->source_side[ixyz]+=cfptr->source_side[ixyz];
					cfbar->source_long[ixyz]+=cfptr->source_long[ixyz];
				}
			}
		}
	}
	for(iq=0;iq<CF::NQ;iq++){
		cfbar->cf_qinv[iq]=cfbar->cf_qinv[iq]/cfbar->norm_qinv[iq];
		cfbar->cf_qout[iq]=cfbar->cf_qout[iq]/cfbar->norm_qout[iq];
		cfbar->cf_qside[iq]=cfbar->cf_qside[iq]/cfbar->norm_qside[iq];
		cfbar->cf_qlong[iq]=cfbar->cf_qlong[iq]/cfbar->norm_qlong[iq];
	}
	printf("nincrement_sum=%lld, cfbar->norm_qinv[0]=%g\n",cfbar->nincrement,cfbar->norm_qinv[0]);
	for(int ictheta=0;ictheta<10;ictheta++){
		for(int iphir=0;iphir<18;iphir++){
			cfbar->ThetaPhiDist[ictheta][iphir]=cfbar->ThetaPhiDist[ictheta][iphir]/double(cfbar->nincrement);
		}
	}
	for(int ixyz=0;ixyz<CF::Nxyz;ixyz++){
		cfbar->source_out[ixyz]+=cfptr->source_out[ixyz]/double(cfbar->nincrement);
		cfbar->source_side[ixyz]+=cfptr->source_side[ixyz]/double(cfbar->nincrement);
		cfbar->source_long[ixyz]+=cfptr->source_long[ixyz]/double(cfbar->nincrement);
	}
}

void CHBT_BES::CalcCoalescenceSpectra(){
	double COAL_SPIN_FACTOR=parmap->getD("COAL_SPIN_FACTOR",0.75);
	double Pt,uperp,D3PoverE,D3PoverEa,D3PoverEb,E,phi,delN,uperp1,uperp0,mu=MASSA*MASSB/(MASSA+MASSB);;
	double Ebar=0.0,Ptbar=0.0,EtPtnorm=0.0,V2=0.0;
	int iuperp,irap,iphi;
	CF *cfptr;
	double cfactor;
	//D3r=4.0*PI*pow(CF::Rcoalescence,3)/3.0; // Coalescence volume in coordinate space
	vector<double> coalspectra,aspectra,bspectra,Rhbt;
	coalspectra.resize(NUPERP);
	aspectra.resize(NUPERP);
	bspectra.resize(NUPERP);
	Rhbt.resize(NUPERP);
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
	for(iuperp=0;iuperp<NUPERP;iuperp++){
		uperp0=iuperp*DELUPERP;
		uperp1=(iuperp+1)*DELUPERP;
		D3PoverEa=DELRAP*NRAP*2.0*PI*MASSA*MASSA*PI*(uperp1*uperp1-uperp0*uperp0);
		D3PoverEb=DELRAP*NRAP*2.0*PI*MASSB*MASSB*PI*(uperp1*uperp1-uperp0*uperp0);
		aspectra[iuperp]=aspectra[iuperp]/(NEVENTS*D3PoverEa);
		bspectra[iuperp]=bspectra[iuperp]/(NEVENTS*D3PoverEb);
		Rhbt[iuperp]=cfactor*aspectra[iuperp]*aspectra[iuperp]/coalspectra[iuperp];
		Rhbt[iuperp]=pow(Rhbt[iuperp],1.0/3.0);
		Pt=(MASSA+MASSB)*(iuperp+0.5)*DELUPERP;
		printf("%7.1f %7.4f\n",Pt,Rhbt[iuperp]);
	}
}

CHBT_Part::CHBT_Part(){
	p.resize(4);
	x.resize(4);
}

void CHBT_Part::Print(){
	printf("ID=%d, mass=%g\n",ID,mass);
	printf("x=(%g,%g,%g,%g), p=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3],p[0],p[1],p[2],p[3]);
}

#endif
