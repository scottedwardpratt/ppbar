#include "hbt_bes.h"

using namespace std;

int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: hbt_ppbar parameters_filename\n");
		exit(1);			
	}
	double Rout,Rside,Rlong;
	
	CHBT_BES hbt(argv[1]);
	printf("Enter Rout,Rside,Rlong: ");
	scanf("%lf %lf %lf",&Rout,&Rside,&Rlong);
	hbt.randy->reset(-time(NULL));
	hbt.CalcCF_Gauss(Rout,Rside,Rlong);
	//hbt.WriteCFs();
	//hbt.WriteThetaPhiDists();
	hbt.cfgauss->Print();
	return 0;
}
