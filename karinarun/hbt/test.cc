#include <stdio.h>      
#include <math.h> 
#include <stdlib.h>
#include <fstream>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"



using namespace std;

int main(int argc,char *argv[]){

	long double x;
	long double r;
	long double arr [10000];
	int ix;
	int delx=20;
	int dndx[200];
	ofstream fout;
	fout.open("test.dat", ios::out);
	
	for(int i=0; i<10000; i++){
		//out of bounds errors?
		r=randy->ran();
		x=-log(r);
		ix=floorl(x/delx);
		if(ix<200){
			dndx[ix]+=1;
			fout << ix << " " << dndx[ix] <<  endl;
		}		
	}
	
    fout << endl;
	
	
	return 0;
}