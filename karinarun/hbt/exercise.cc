#include <stdio.h>      
#include <math.h> 
#include <stdlib.h>
#include <fstream>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"
#include "randy.h"



using namespace std;

int main(int argc,char *argv[]){

	double x;
	double r;
	CRandy *randy=new CRandy(-12345);
	double arr [10000];
	int ix;
	int delx=20;
	int dndx[200];
	ofstream fout;
	double T=300.0;
	fout.open("exercise.dat", ios::out);
	
	for(int i=0; i<1000000; i++){
		r=randy->ran();
//      printf("r: %g\n", r);
		x=-T*log(r);
		ix=floorl(x/delx);		
			if(ix<200){
			dndx[ix]+=1;
		}		
	}
	for(int i=0; i<200; i++){
		fout << i << " " << dndx[i] <<  endl;
		
	}
	
    fout << endl;
	
	
	return 0;
}