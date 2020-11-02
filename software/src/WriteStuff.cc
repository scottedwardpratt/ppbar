#ifndef __HBT_BES_WRITESTUFF_CC__
#define __HBT_BES_WRITESTUFF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"


void CHBT_BES::WriteCFs(){
	char filename[150];
	int irap,iphi,ipt;
	double DELPHI=360.0/double(NPHI);
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			for(ipt=0;ipt<NPT;ipt++){
				sprintf(filename,"results/CFs/rap%g_phi%g_pt%g.txt",(irap+0.5)*DELRAP,(iphi+0.5)*DELPHI,(ipt+0.5)*DELPT);
				CFArray[irap][iphi][ipt]->WriteCFs(string(filename));
			}
		}
	}
	cfbar->WriteCFs("results/CFs/average.txt");
}

void CHBT_BES::WriteThetaPhiDists(){
	char filename[150];
	int irap,iphi,ipt;
	double DELPHI=360.0/double(NPHI);
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			for(ipt=0;ipt<NPT;ipt++){
				sprintf(filename,"results/ThetaPhiDists/rap%g_phi%g_pt%g.txt",(irap+0.5)*DELRAP,(iphi+0.5)*DELPHI,(ipt+0.5)*DELPT);
				CFArray[irap][iphi][ipt]->WriteThetaPhiDist(string(filename));
			}
		}
	}
	cfbar->WriteThetaPhiDist("results/ThetaPhiDists/average.txt");
}

#endif