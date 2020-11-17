#ifndef __HBT_BES_CALC_CF_CC__
#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

void CHBT_BES::CalcCF_MC(){
	CHBT_PartMap::iterator ita,itb;
	int irap,iphi,iuperp;
	vector<double> x(4);
	CHBT_Part *partaa,*partbb;
	CF *cf;
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			printf("irap=%d, iphi=%d, size=%d\n",irap,iphi,int(partmap[irap][iphi].size()));
			for(ita=partmap[irap][iphi].begin();ita!=partmap[irap][iphi].end();++ita){
				partaa=ita->second;
				if(partaa->ID==IDA || partaa->ID==IDB){
					itb=ita; itb++;
					partbb=itb->second;
					while(itb!=partmap[irap][iphi].end() && fabs(partaa->uperp-partbb->uperp)<UPERPTEST){
						partbb=itb->second;
						if((partaa->ID==IDA && partbb->ID==IDB) || (partaa->ID==IDB && partbb->ID==IDA)){
							cf=GetCF(partaa,partbb);
							if(cf!=NULL){
								if(partaa->ID==IDA){
									cf->Increment(partaa,partbb);
								}
								else
									cf->Increment(partbb,partaa);
							}
						}
						itb++;
					}
				}
			}
		}
	}
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			for(iuperp=0;iuperp<NUPERP;iuperp++){
				CFArray[irap][iphi][iuperp]->Normalize();
			}
		}
	}
}

void CHBT_BES::CalcCF_Gauss(double Rout,double Rside,double Rlong){
	int imc;
	vector<double> x;
	x.resize(4);
	x[0]=0.0;
	cfgauss->nincrement=0;
	for(imc=0;imc<NMC;imc++){
		x[1]=Rout*randy->ran_gauss();
		x[2]=Rside*randy->ran_gauss();
		x[3]=Rlong*randy->ran_gauss();
		cfgauss->Increment(x);
		if((10*(imc+1))%NMC==0){
			printf("finished %g percent\n",100.0*(imc+1)/double(NMC));
		}
	}
	cfgauss->Normalize();
}


#endif