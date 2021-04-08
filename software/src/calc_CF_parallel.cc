//#ifndef __HBT_BES_CALC_CF_CC__
//#define __HBT_BES_CALC_CF_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"
#include "omp.h"
void CHBT_BES::CalcCF(){
	CHBT_PartMap::iterator ita,itb;
	int irap,iphi,iuperp;
	vector<double> x(4);
	CHBT_Part *partaa,*partbb;
	CF *cf;
	long long int ntot=0,ntotprime=0;
	for(irap=0;irap<NRAP;irap++){
		for(iphi=0;iphi<NPHI;iphi++){
			ntot+=partmap[irap][iphi].size();
			printf("irap=%d, iphi=%d, size=%d, ntot=%lld\n",irap,iphi,int(partmap[irap][iphi].size()),ntot);
			for(ita=partmap[irap][iphi].begin();ita!=partmap[irap][iphi].end();++ita){
				partaa=ita->second;
				ntotprime+=1;
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

void CHBT_BES::CalcCF_Gauss(double Rout,double Rside,double Rlong, double xoff, double yoff, double zoff){
	int imc;
	vector<double> x;
	double root2=sqrt(2.0);
	x.resize(4);
	x[0]=0.0;
	int j=1;
	cfgauss->nincrement=0;
	double x2,y,z;
	//#pragma omp parallel for
	for(imc = 0;imc<NMC;imc++){
		x2=root2*Rout*randy->ran_gauss();
		y=root2*Rside*randy->ran_gauss();
		z=root2*Rlong*randy->ran_gauss();
		
		x[1]=x2+xoff;
		x[2]=y+yoff;
		x[3]=z+zoff;
		cfgauss->Increment(x);
		
		x[1]=-x2+xoff;
		x[2]=-y+yoff;
		x[3]=-z+zoff;
		cfgauss->Increment(x);
			
	}
	
				
	
	/*
	if((10*(imc+1))%NMC==0){
	printf("finished %g percent\n",100.0*(imc+1)/double(NMC));
	}
	*/
	
cfgauss->Normalize();
}


//#endif