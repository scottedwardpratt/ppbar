#include "hbt_bes.h"

using namespace std;

int main(){
	CHBT_BES hbt("parameters.txt");
	hbt.randy->reset(-time(NULL));
	hbt.ReadPR();
	hbt.CalcCF_MC();
	//hbt.CFArray[0][0][10]->Print();
	hbt.AverageCF();
	hbt.cfbar->Print();
	hbt.WriteCFs();
	hbt.WriteThetaPhiDists();
	return 0;
}
