#include "hbt_bes.h"
#include "time.h"
using namespace std;

int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: hbt_ppbar parameters_filename\n");
		exit(1);			
	}
	double Rout,Rside,Rlong, xoff, yoff, zoff;
	CHBT_BES hbt(argv[1]);
	hbt.randy->reset(-time(NULL));
	CF cf_smash;
	string filename;
	time_t start, end; 
	
	filename=hbt.RESULTS_DIR+"/CFs/average.txt";
	cf_smash.ReadCFs(filename);

	printf("Enter Rout,Rside,Rlong, xoff, yoff, zoff: ");
	scanf("%lf %lf %lf %lf %lf %lf",&Rout,&Rside,&Rlong,&xoff,&yoff,&zoff);
	hbt.RESULTS_DIR="results_gauss";
	start=time(NULL);
	hbt.CalcCF_Gauss(Rout,Rside,Rlong,xoff,yoff,zoff);
    end=time(NULL);
	
	double chisquare=hbt.GetChiSquare(hbt.cfgauss,&cf_smash);
	printf("chisquare=%g\n",chisquare);
	//printf("Time elapsed:%.2f\n", difftime(end,start));
	
	hbt.cfgauss->Print();
	return 0;
}
