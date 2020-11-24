#include "hbt_bes.h"

using namespace std;

int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: hbt_ppbar parameters_filename\n");
		exit(1);			
	}
	double Rout,Rside,Rlong;
	CHBT_BES hbt(argv[1]);
	hbt.randy->reset(-time(NULL));
	CF cf_smash;
	string filename;
	
	filename=hbt.RESULTS_DIR+"/CFs/average.txt";
	cf_smash.ReadCF(filename);

	printf("Enter Rout,Rside,Rlong: ");
	scanf("%lf %lf %lf",&Rout,&Rside,&Rlong);
	hbt.RESULTS_DIR="results_gauss";
	hbt.CalcCF_Gauss(Rout,Rside,Rlong);
	
	double chisquare=hbt.GetChiSquare(hbt.cfgauss,&cf_smash);
	printf("chisquare=%g\n",chisquare);
	hbt.cfgauss->Print();
	return 0;
}
