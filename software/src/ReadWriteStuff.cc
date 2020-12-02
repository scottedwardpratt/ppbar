#ifndef __HBT_BES_READWRITESTUFF_CC__
#define __HBT_BES_READWRITESTUFF_CC__
#include <cstdlib>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

void CHBT_BES::ReadPR(){
	int iskip,count;
	int ipart,nparts,ID,charge,smashID,ievent,irun;
	char dummy[200],dumbo1[20],dumbo2[20],dumbo3[20];
	double ptbar=0.0,ptnorm=0.0,pt,rap;
	vector<double> p,x;
	p.resize(4);
	x.resize(4);
	double mass;
	NEVENTS=0;
	FILE *oscarfile;
	for(irun=0;irun<INPUT_OSCAR_NRUNS;irun++){
		INPUT_OSCAR_FILENAME=INPUT_OSCAR_BASE_DIRECTORY+"run"+to_string(irun)+"/"+"0/particle_lists.oscar";
		oscarfile=fopen(INPUT_OSCAR_FILENAME.c_str(),"r");
		printf("opening %s\n",INPUT_OSCAR_FILENAME.c_str());
		for(iskip=0;iskip<3;iskip++){
			fgets(dummy,200,oscarfile);
		}
		while(!feof(oscarfile) && NEVENTS<NEVENTS_MAX){
			fscanf(oscarfile,"%s %s %d %s %d",dumbo1,dumbo2,&ievent,dumbo3,&nparts);
			//printf("nparts=%d\n",nparts);
			fgets(dummy,160,oscarfile);
			if(!feof(oscarfile)){
				NEVENTS+=1;
				if(NEVENTS==NEVENTS_MAX)
					irun=INPUT_OSCAR_NRUNS;
				count=0;
				for(ipart=0;ipart<nparts;ipart++){
					fscanf(oscarfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
					&x[0],&x[1],&x[2],&x[3],&mass,&p[0],&p[1],&p[2],&p[3],&ID,&smashID,&charge);
					fgets(dummy,160,oscarfile);
					//printf("check ipart=%d, ID=%d\n",ipart,ID);
					if(ID==IDA){
						pt=1000.0*sqrt(p[1]*p[1]+p[2]*p[2]);
						rap=asinh(p[3]*1000.0/sqrt(1000.0*mass*1000.0*mass+pt*pt));
						AddPart(ID,p,x,mass);
						if(fabs(rap)<0.5*NRAP*DELRAP){
							ptnorm+=1.0;
							ptbar+=pt;
						}
					}
					else if(IDB!=IDA && ID==IDB){
						AddPart(ID,p,x,mass);
					}
					count+=1;
				}
				fgets(dummy,200,oscarfile);
			}
		}
		fclose(oscarfile);
	}
	printf("for particle A, <pt>=%g, N_A=%g NEVENTS=%d\n",ptbar/ptnorm,ptnorm,NEVENTS);
}

void CHBT_BES::WriteCFs(){
	char filename[150];
	int irap,iphi,iuperp;
	char command[120];
	sprintf(command,"mkdir -p %s/CFs",RESULTS_DIR.c_str());
	system(command);
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			for(iuperp=0;iuperp<NUPERP;iuperp++){
				sprintf(filename,"%s/CFs/rap%g_phi%g_Pt%g.txt",RESULTS_DIR.c_str(),
				-0.5*NRAP*DELRAP+(irap+0.5)*DELRAP,-0.5*NPHI*DELPHI*180.0/PI+(iphi+0.5)*DELPHI*180.0/PI,(iuperp+0.5)*DELPT);
				CFArray[irap][iphi][iuperp]->WriteCFs(string(filename));
			}
		}
	}
	sprintf(filename,"%s/CFs/average.txt",RESULTS_DIR.c_str());
	cfbar->WriteCFs(string(filename));
}

void CHBT_BES::ReadCFs(){
	char filename[150];
	int irap,iphi,iuperp;
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			for(iuperp=0;iuperp<NUPERP;iuperp++){
				sprintf(filename,"%s/CFs/rap%g_phi%g_Pt%g.txt",RESULTS_DIR.c_str(),
				-0.5*NRAP*DELRAP+(irap+0.5)*DELRAP,-0.5*NPHI*DELPHI*180.0/PI+(iphi+0.5)*DELPHI*180.0/PI,(iuperp+0.5)*DELPT);
				CFArray[irap][iphi][iuperp]->ReadCFs(string(filename));
			}
		}
	}
}

#endif