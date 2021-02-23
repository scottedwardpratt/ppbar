#include "hbt_bes.h"
#include <omp.h>
//#include "mpi.h"
#include "time.h"

using namespace std;

int main(int argc,char *argv[]){
	const int NTHREADS=24;
	//int rank, numtasks;
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
	CHBT_BES hbt_fake(argv[1]);
	hbt_fake.randy->reset(-time(NULL));
	//CF cf_smash;
	//	if(rank==0){
	string filename;
	filename=hbt_fake.RESULTS_DIR+"/CFs/average.txt";
	hbt_fake.CalcCF_Gauss(3.0,4.0,5.0);
	//cf_smash.ReadCFs(filename);
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
	double best[3]={3.5, 3.5, 5.0};
	double diff=1.0,chisquare_best=5600000.0;
	int location;
	for(int k=0; k<100; ++k){
		printf("k=%d, best R=(%g,%g,%g), diff=%g\n",k,best[0],best[1],best[2],diff);
		
		
#pragma omp parallel for
		for(int i=0; i<NTHREADS; ++i){
			//my innovating
			if(i==0){
				for(int j=0; j<3; j++)
					arr[i][j]=best[j];
			}
			else{
				for(int j=0; j<3; j++){
					arr[i][j]=best[j]+diff*(1.0-2.0*hbt_gauss[i]->randy->ran());
				}
			}
			hbt_gauss[i]->cfgauss->Reset();
			hbt_gauss[i]->CalcCF_Gauss(arr[i][0],arr[i][1],arr[i][2]);
			chisquare[i]=hbt_gauss[i]->GetChiSquare(hbt_gauss[i]->cfgauss,hbt_fake.cfgauss);	
		}
		
		location=-1;
		bool success=false;
		chisquare_best=100000000.0;
		for(int i=0;i<NTHREADS;i++){
			if(chisquare[i]<chisquare_best){
				if(i!=0)
					success=true;
				chisquare_best=chisquare[i];
				location=i;
				for(int j=0; j<3; j++){
					best[j]= arr[location][j];
				}
			}
		}
		if(success){
			printf("success, diff=%g\n",diff);
			printf("chisquare=%g\n",chisquare_best);
			printf("Rout:%g\n", best[0]);
			printf("Rside:%g\n", best[1]);
			printf("Rlong:%g\n", best[2]);
			hbt_gauss[location]->cfgauss->Print();
		}
		if(!success)
			diff*=0.75;
		//printf("finished thread %d\n",i);
	}
	printf("Target CF:\n");
	hbt_fake.cfgauss->Print();

	//	if(rank==0){
	end=time(NULL);
	printf("Time elapsed:%g\n", difftime(end,start));
	
	//	}
	return 0;
	//	MPI_Finalize();
}
