#include "hbt_bes.h"

using namespace std;

int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: hbt_ppbar parameters_filename\n");
		exit(1);			
	}
	double ptmin=1000.0,ptmax=2000.0,rapmin=-0.6,rapmax=0.6,phimin_deg=-180,phimax_deg=180;
	double Rout,Rside,Rlong;
	CHBT_BES hbt_smash(argv[1]);
	hbt_smash.randy->reset(-time(NULL));
	hbt_smash.ReadCFs();
	hbt_smash.AverageCF(hbt_smash.cfbar,rapmin,rapmax,ptmin,ptmax,phimin_deg,phimax_deg);
	hbt_smash.cfbar->Print();

	printf("Enter Rout,Rside,Rlong: ");
	scanf("%lf %lf %lf",&Rout,&Rside,&Rlong);
	hbt_smash.RESULTS_DIR="results_gauss";
	hbt_smash.CalcCF_Gauss(Rout,Rside,Rlong);
	
	double chisquare=hbt_smash.GetChiSquare(hbt_smash.cfgauss,hbt_smash.cfbar);
	printf("chisquare=%g\n",chisquare);
	hbt_smash.cfgauss->Print();
	
	return 0;
}
