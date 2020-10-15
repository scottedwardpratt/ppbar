#ifndef __HBT_BES_CC__
#define __HBT_BES_CC__
#include <cstring>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

CHBT_BES::CHBT_BES(string parsfilename){
	parmap=new CparameterMap();
	parmap->ReadParsFromFile(parsfilename);
	TAU_COMPARE=parmap->getD("TAU_COMPARE",12.0);
	INPUT_OSCAR_FILENAME=parmap->getS("INPUT_OSCAR_FILENAME","data/crap.oscar");
	Nqinv=parmap->getI("Nqinv",100);
	DELqinv=parmap->getD("DELqinv",1.0);
	qinvMAX=Nqinv*DELqinv;
	IDA=parmap->getI("IDA",-2212);
	IDB=parmap->getI("IDB",-2212);
	NEVENTS_MAX=parmap->getI("NEVENTS_MAX",10);
	NMC=parmap->getI("NMC",100000);
	INPUT_OSCAR_BASE_DIRECTORY=parmap->getS("INPUT_OSCAR_BASE_DIRECTORY","crap");
	INPUT_OSCAR_NRUNS=parmap->getI("INPUT_OSCAR_NRUNS",24);
	RANSEED=parmap->getD("RANSEED",-12345);
	QINVTEST=parmap->getD("QINVTEST",50.0);
	randy=new CRandy(RANSEED);
	
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
	CFqinv.resize(Nqinv);
	
}

void CHBT_BES::ReadPR(){
	double YMAX=0.5;
	int ipart,nparts,ID,nevents,charge,smashID,irun;
	char dummy[200],dumbo1[20],dumbo2[20],dumbo3[20];
	vector<double> p,x;
	double phi,yy,sphi,cphi,rap,gamma,gammav,pt,z,taucalc=25.0;
	p.resize(4);
	x.resize(4);
	double mass;
	FILE *oscarfile;
	for(irun=0;irun<INPUT_OSCAR_NRUNS;irun++){
		INPUT_OSCAR_FILENAME=INPUT_OSCAR_BASE_DIRECTORY+"run"+to_string(irun)+"/"+"0/particle_lists.oscar";
		oscarfile=fopen(INPUT_OSCAR_FILENAME.c_str(),"r");
		printf("opening %s\n",INPUT_OSCAR_FILENAME.c_str());
		for(int idummy=0;idummy<3;idummy++){
			fgets(dummy,200,oscarfile);
		}
		do{
			fscanf(oscarfile,"%s %s %d %s %d",dumbo1,dumbo2,&nevents,dumbo3,&nparts);
			printf("nparts=%d\n",nparts);
			fgets(dummy,160,oscarfile);
			if(!feof(oscarfile)){
				for(ipart=0;ipart<nparts;ipart++){
					fscanf(oscarfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
					&x[0],&x[1],&x[2],&x[3],&mass,&p[0],&p[1],&p[2],&p[3],&ID,&smashID,&charge);
					fgets(dummy,160,oscarfile);
					//printf("check ipart=%d, ID=%d\n",ipart,ID);
					for(int alpha=0;alpha<4;alpha++){
						p[alpha]*=1000.0;
						mass*=1000.0;
					}
					rap=atanh(p[3]/p[0]);
					if(fabs(rap)<YMAX){
						phi=atan2(p[2],p[1]);
						pt=sqrt(p[1]*p[1]+p[2]*p[2]);
						cphi=cos(phi);
						sphi=sin(phi);
						//r=sqrt(x[1]*x[1]+x[2]*x[2]);
						yy=x[2];
						x[2]=x[2]*cphi-x[1]*sphi;
						x[1]=x[1]*cphi+yy*sphi;
						p[2]=0.0;
						p[1]=pt;
						gamma=cosh(rap);
						gammav=sinh(rap);
						p[0]=p[0]*gamma-gammav*p[3];
						p[3]=0.0;
						z=p[3];
						x[3]=gamma*x[3]-gammav*x[0];
						x[0]=gamma*x[0]-gammav*z;
						x[3]=x[3]-(p[3]/p[0])*(x[0]-taucalc);
						x[0]=taucalc;
						if(ID==IDA){
							AddPart(parta,ID,p,x);
						}
						if(IDB!=IDA && ID==IDB){
							AddPart(partb,ID,p,x);
						}
					}
				}
				fgets(dummy,160,oscarfile);
			}
		}while(!feof(oscarfile) && nevents<NEVENTS_MAX);
		fclose(oscarfile);
	}
}

void CHBT_BES::AddPart(vector<CHBT_Part *> &part,int &IDread,vector<double> &pread,vector<double> &xread){
	CHBT_Part *newpart;
	newpart=new CHBT_Part();
	part.push_back(newpart);
	newpart->ID=IDread;
	for(int alpha=0;alpha<4;alpha++){
		newpart->p[alpha]=pread[alpha];
	}
	CalcXBjPt(newpart,xread);
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

void CHBT_BES::CalcXR(CHBT_Part *partaa,CHBT_Part *partbb,vector<double> &x,double &r){
	double dummy,gamma,gammav=0.5*((partaa->p[1]/partaa->mass)+(partbb->p[1]/partbb->mass));
	gamma=sqrt(1.0+gammav*gammav);
	for(int alpha=0;alpha<4;alpha++)
		x[alpha]=partaa->xbj[alpha]-partbb->xbj[alpha];
	dummy=x[0];
	x[0]=gamma*dummy-gammav*x[1];
	x[1]=gamma*x[1]-gammav*dummy;
	r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
}


void CHBT_BES::CalcXBjPt(CHBT_Part *part,vector<double> &x){
	double gamma,gammav,mperp,phi,cphi,sphi;
	mperp=sqrt(part->p[0]*part->p[0]-part->p[3]*part->p[3]);
	part->pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	part->mass=sqrt(part->p[0]*part->p[0]-part->pt*part->pt-part->p[3]*part->p[3]);
	gammav=part->p[3]/mperp;
	gamma=sqrt(1.0+gammav*gammav);
	part->xbj[0]=gamma*x[0]-gammav*x[3];
	part->xbj[3]=gamma*x[3]-gammav*x[0];
	phi=atan2(part->p[2],part->p[1]);
	cphi=cos(phi);
	sphi=sin(phi);
	part->xbj[1]=x[1]*cphi+x[2]*sphi;
	part->xbj[2]=x[2]*cphi-x[1]*sphi;
	part->xbj[1]=part->xbj[1]-(part->xbj[0]-TAU_COMPARE)*part->pt/mperp;
	part->xbj[0]=TAU_COMPARE;	
}

CHBT_Part::CHBT_Part(){
	p.resize(4);
	xbj.resize(4);
}
void CHBT_Part::Print(){
	printf("ID=%d, mass=%g\n",ID,mass);
	printf("xbj=(%g,%g,%g,%g), p=(%g,%g,%g,%g)\n",xbj[0],xbj[1],xbj[2],xbj[3],p[0],p[1],p[2],p[3]);
}
#endif
