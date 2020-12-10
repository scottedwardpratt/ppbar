#include "hbt_bes.h"
#include <omp.h>
//#include "mpi.h"
#include "time.h"

using namespace std;

int main(int argc,char *argv[]){
	const int NTHREADS=24;
	int rank, numtasks;
	//MPI initialization
	/*	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	*/	
	omp_set_num_threads(NTHREADS);
	
	if(argc!=2){
		printf("Usage: hbt_ppbar parameters_filename\n");
		exit(1);			
	}
	double chisquare[NTHREADS];
	CHBT_BES *hbt_gauss[NTHREADS];
	CHBT_BES hbt_smash(argv[1]);
	hbt_smash.randy->reset(-time(NULL));
	CF cf_smash;
	//	if(rank==0){
	string filename;
	filename=hbt_smash.RESULTS_DIR+"/CFs/average.txt";
	cf_smash.ReadCFs(filename);
	//	}
	//	MPI_Scatter(&cf_smash, 250, MPI_LONG_DOUBLE, &cf_smash, 250, MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);
	time_t start,end;
	start=time(NULL);
	for(int i=0; i<NTHREADS; ++i){
		char dirname[100];
		char command[120];
		sprintf(dirname,"results_gauss%d",i);
		sprintf(command,"mkdir -p %s\n",dirname);
		system(command);
		hbt_gauss[i]=new CHBT_BES(argv[1]);
		hbt_gauss[i]->randy->reset(i);
		hbt_gauss[i]->RESULTS_DIR=dirname;
	}
	

	double arr[NTHREADS][3];
	double best[3]={5.0, 5.0, 5.0};
	double diff=3.0;
	int location;
	for(int k=0; k<10; ++k){
		
		
#pragma omp parallel for
		for(int i=0; i<NTHREADS; ++i){
			//my innovating
			for(int j=0; j<3; j++){
				//	best[j]=5.0;
				arr[i][j]=best[j]+diff*(1.0-2.0*hbt_gauss[i]->randy->ran());
			}
			double chisquare1=560000.0;
		
			hbt_gauss[i]->CalcCF_Gauss(arr[i][0],arr[i][1],arr[i][2]);
			chisquare[i]=hbt_gauss[i]->GetChiSquare(hbt_gauss[i]->cfgauss,&cf_smash);	
			location=-1;
			if(chisquare[i]<chisquare1){
				location=i;
			}
			printf("diff value:%g\n", diff);
			printf("finished thread %d\n",i);
		}
								
				
		//end of innovating
		//part i need
		/*		for(int j=0; j<3; j++){
		arr[i][j]=2.0+8.0*hbt_gauss[i]->randy->ran();
		}
		hbt_gauss[i]->CalcCF_Gauss(arr[i][0],arr[i][1],arr[i][2]);
		chisquare[i]=hbt_gauss[i]->GetChiSquare(hbt_gauss[i]->cfgauss,&cf_smash);	
		*/		//end of part i need

		//		double minchisquare[12];
		//	MPI_Reduce(&chisquare, &minchisquare, 1000, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);	
		//   double chisquare=hbt.GetChiSquare(hbt.cfgauss,&cf_smash);
		//finds minimum value in minchisquar which returns an array of size 12
/*		for(int i=1; i<NTHREADS; i++){
			if(chisquare[i]<chisquare[location]){
				location=i;
			}
		}
*/		for(int j=0; j<3; j++){
			best[j]= arr[location][j];
		}
		diff=0.5*diff;
		printf("chisquare=%g\n",chisquare[location]);
		printf("Rout:%g\n", arr[location][0]);
		printf("Rside:%g\n", arr[location][1]);	
		printf("Rlong:%g\n", arr[location][2]);
	}

	//	if(rank==0){
	end=time(NULL);
	printf("Time elapsed:%g\n", difftime(end,start));
	
	for(int i=0;i<NTHREADS;i++){
		hbt_gauss[i]->cfgauss->Print();
	}
	//	}
	return 0;
	//	MPI_Finalize();
}
