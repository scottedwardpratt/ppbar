#ifndef __HBT_BES_WRITESTUFF_CC__
#define __HBT_BES_WRITESTUFF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"


void CHBT_BES::WriteCFs(){
	char filename[150];
	int irap,iphi,iuperp;
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			for(iuperp=0;iuperp<NUPERP;iuperp++){
				sprintf(filename,"%s/CFs/rap%g_phi%g_uperp%g.txt",RESULTS_DIR.c_str(),
				(irap+0.5)*DELRAP,(iphi+0.5)*DELPHI*180.0/PI,(iuperp+0.5)*DELUPERP);
				CFArray[irap][iphi][iuperp]->WriteCFs(string(filename));
			}
		}
	}
	sprintf(filename,"%s/CFs/average.txt",RESULTS_DIR.c_str());
	cfbar->WriteCFs(string(filename));
}

void CHBT_BES::WriteThetaPhiDists(){
	char filename[150];
	int irap,iphi,iuperp;
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			for(iuperp=0;iuperp<NUPERP;iuperp++){
				sprintf(filename,"%s/ThetaPhiDists/rap%g_phi%g_uperp%g.txt",RESULTS_DIR.c_str(),(irap+0.5)*DELRAP,(iphi+0.5)*DELPHI,(iuperp+0.5)*DELUPERP);
				CFArray[irap][iphi][iuperp]->WriteThetaPhiDist(string(filename));
			}
		}
	}
	sprintf(filename,"%s/ThetaPhiDists/average.txt",RESULTS_DIR.c_str());
	cfbar->WriteThetaPhiDist(string(filename));
}

#endif