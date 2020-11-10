#ifndef __HBT_BES_CC__
#define __HBT_BES_CC__
#include <cstring>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

CHBT_BES::CHBT_BES(string parsfilename){
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
	RESULTS_DIR=parmap->getS("RESULTS_DIR","results");
	NPHI=parmap->getI("NPHI",16);
	NRAP=parmap->getI("NRAP",3);
	DELRAP=parmap->getD("DELRAP",0.2);
	NPT=parmap->getI("NPT",20);
	DELPT=parmap->getD("DELPT",100.0);
	YMAX=1.0;
	IDA=parmap->getI("IDA",211);
	IDB=parmap->getI("IDB",211);
	NEVENTS_MAX=parmap->getI("NEVENTS_MAX",10);
	NMC=parmap->getI("NMC",100000);
	TAU_COMPARE=25.0;
	INPUT_OSCAR_BASE_DIRECTORY=GITHOME_MSU+"/"+parmap->getS("INPUT_OSCAR_BASE_DIRECTORY","crap");
	INPUT_OSCAR_NRUNS=parmap->getI("OSCAR_NRUNS",24);
	//RANSEED=parmap->getD("RANSEED",-12345);
	RANSEED=-time(NULL);
	QINVTEST=parmap->getD("QINVTEST",50.0);
	randy=new CRandy(RANSEED);
	CF::randy=randy;
	CF::hbt=this;
	
	if((IDA==2212 && IDB==2212) || (IDA==-2212 && IDB==-2212)){
		wf=new CWaveFunction_pp_schrod(parsfilename);
	}
	else if((IDA==211 && IDB==211) || (IDA==-211 && IDB==-211)){
		wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);
	}
	else if((IDA==211 && IDB==-211) || (IDA==-211 && IDB==211)){
		wf=new CWaveFunction_pipluspiminus_sqwell(parsfilename);
	}
	else if((IDA==321 && IDB==321) || (IDA==-321 && IDB==-321)){
		wf=new CWaveFunction_generic(parsfilename,1,493.677,493.677,1.0);
	}
	else if((IDA==321 && IDB==-321) || (IDA==-321 && IDB==321)){
		wf=new CWaveFunction_generic(parsfilename,-1,493.677,493.677,1.0);
	}
	else{
		printf("fatal: IDs not recognized, IDA=%d, IDB=%d\n",IDA,IDB);
		exit(1);
	}
	
	if(!GAUSS){
		CFArray.resize(NRAP);
		for(int irap=0;irap<NRAP;irap++){
			CFArray[irap].resize(NPHI);
			for(int iphi=0;iphi<NPHI;iphi++){
				CFArray[irap][iphi].resize(NPT);
				for(int ipt=0;ipt<NPT;ipt++){
					CFArray[irap][iphi][ipt]=new CF();
					CFArray[irap][iphi][ipt]->Reset();
					CFArray[irap][iphi][ipt]->wf=wf;
				}
			}
		}
		cfbar=new CF();
		cfbar->Reset();
		cfbar->wf=wf;
	}
	else{
		cfgauss=new CF();
		cfgauss->nincrement=0;
		cfgauss->Reset();
		cfgauss->wf=wf;
	}
}

void CHBT_BES::ReadPR(){
	int iskip,count;
	int ipart,nparts,ID,nevents,charge,smashID,irun;
	char dummy[200],dumbo1[20],dumbo2[20],dumbo3[20];
	double rap;
	vector<double> p,x;
	p.resize(4);
	x.resize(4);
	double mass;
	FILE *oscarfile;
	for(irun=0;irun<INPUT_OSCAR_NRUNS;irun++){
		INPUT_OSCAR_FILENAME=INPUT_OSCAR_BASE_DIRECTORY+"run"+to_string(irun)+"/"+"0/particle_lists.oscar";
		oscarfile=fopen(INPUT_OSCAR_FILENAME.c_str(),"r");
		printf("opening %s\n",INPUT_OSCAR_FILENAME.c_str());
		for(iskip=0;iskip<3;iskip++){
			fgets(dummy,200,oscarfile);
		}
		do{
			fscanf(oscarfile,"%s %s %d %s %d",dumbo1,dumbo2,&nevents,dumbo3,&nparts);
			//printf("nparts=%d\n",nparts);
			fgets(dummy,160,oscarfile);
			if(!feof(oscarfile)){
				count=0;
				for(ipart=0;ipart<nparts;ipart++){
					fscanf(oscarfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
					&x[0],&x[1],&x[2],&x[3],&mass,&p[0],&p[1],&p[2],&p[3],&ID,&smashID,&charge);
					fgets(dummy,160,oscarfile);
					//printf("check ipart=%d, ID=%d\n",ipart,ID);
					if(ID==IDA){
						rap=atanh(p[3]/p[0]);
						if(fabs(rap)<YMAX)
							AddPart(parta,ID,p,x);
					}
					else if(IDB!=IDA && ID==IDB){
						rap=atanh(p[3]/p[0]);
						if(fabs(rap)<YMAX)
							AddPart(parta,ID,p,x);
					}
					count+=1;
				}
				fgets(dummy,200,oscarfile);
			}
		}while(!feof(oscarfile) && nevents<NEVENTS_MAX);
		fclose(oscarfile);
	}
}

void CHBT_BES::AddPart(vector<CHBT_Part *> &part,int &IDread,vector<double> &p,vector<double> &x){
	CHBT_Part *newpart;
	double rap0,phi0,pt,y,z,cphi,sphi,gamma,gammav;
	newpart=new CHBT_Part();
	part.push_back(newpart);
	newpart->ID=IDread;
	for(int alpha=0;alpha<4;alpha++){
		p[alpha]*=1000.0;
	}
	rap0=atanh(p[3]/p[0]);
	
	phi0=atan2(p[2],p[1]);
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	cphi=cos(phi0);
	sphi=sin(phi0);
	//r=sqrt(x[1]*x[1]+x[2]*x[2]);
	y=x[2];
	x[2]=x[2]*cphi-x[1]*sphi;
	x[1]=x[1]*cphi+y*sphi;
	p[2]=0.0;
	p[1]=pt;
	gamma=cosh(rap0);
	gammav=sinh(rap0);
	p[0]=p[0]*gamma-gammav*p[3];
	p[3]=0.0;
	z=x[3];
	x[3]=gamma*x[3]-gammav*x[0];
	x[0]=gamma*x[0]-gammav*z;
	x[1]=x[1]-(p[1]/p[0])*(x[0]-TAU_COMPARE);
	x[0]=TAU_COMPARE;
	for(int alpha=0;alpha<4;alpha++){
		newpart->p[alpha]=p[alpha];
		newpart->x[alpha]=x[alpha];
	}
	newpart->mass=sqrt(p[0]*p[0]-pt*pt-p[3]*p[3]);
	newpart->pt=pt;
	newpart->rap0=rap0;
	newpart->phi0=phi0;
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
	double phi,pt,rap,pxa,pya,pxb,pyb,PX,PY;
	int irap,iphi,ipt;
	rap=0.5*(partaa->rap0+partbb->rap0);
	irap=fabs(rap)/DELRAP;
	if(irap<NRAP){
		pxa=partaa->pt*cos(partaa->phi0);
		pya=partaa->pt*sin(partaa->phi0);
		pxb=partbb->pt*cos(partbb->phi0);
		pyb=partbb->pt*sin(partbb->phi0);
		PX=pxa+pxb;
		PY=pya+pyb;
		phi=atan2(PY,PX);
		pt=0.5*sqrt(PX*PX+PY*PY);
		ipt=pt/DELPT;
		if(ipt<NPT){
			if(phi<0.0)
				phi+=2.0*PI;
			iphi=lrint(floor(phi*NPHI/(2.0*PI)));
			if(iphi<0 || iphi>=NPHI || ipt<0 || ipt>=NPT || irap<0 || irap>=NRAP){
				printf("oopsie\n");
				exit(1);
			}
			return CFArray[irap][iphi][ipt];
		}
		else
			return NULL;
	}
	else
		return NULL;
}

void CHBT_BES::AverageCF(){
	CF *cfptr;
	int iq;
	//long long int cfbar->norm_thetaphidist=0;
	cfbar->Reset();
	for(int irap=0;irap<NRAP;irap++){
		for(int iphi=0;iphi<NPHI;iphi++){
			for(int ipt=0;ipt<NPT;ipt++){
				cfptr=CFArray[irap][iphi][ipt];
				cfbar->nincrement+=cfptr->nincrement;
				if(cfptr->norm_qinv[0]>0){
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
							printf("irap=%d, iphi=%d, ipt=%d\n",irap,iphi,ipt);
							CFArray[irap][iphi][ipt]->Print();
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
		cfbar->cf_qinv[iq]=cfbar->cf_qinv[iq]/double(cfbar->norm_qinv[iq]);
		cfbar->cf_qout[iq]=cfbar->cf_qout[iq]/double(cfbar->norm_qout[iq]);
		cfbar->cf_qside[iq]=cfbar->cf_qside[iq]/double(cfbar->norm_qside[iq]);
		cfbar->cf_qlong[iq]=cfbar->cf_qlong[iq]/double(cfbar->norm_qlong[iq]);
	}
	printf("nincrement_sum=%lld, cfbar->norm_qinv[0]=%g\n",cfbar->nincrement,double(cfbar->norm_qinv[0]));
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

CHBT_Part::CHBT_Part(){
	p.resize(4);
	x.resize(4);
}

void CHBT_Part::Print(){
	printf("ID=%d, mass=%g\n",ID,mass);
	printf("x=(%g,%g,%g,%g), p=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3],p[0],p[1],p[2],p[3]);
}

#endif
