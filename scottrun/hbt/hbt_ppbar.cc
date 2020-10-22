#include "hbt_bes.h"

using namespace std;

int main(){
	CHBT_BES hbt("parameters.txt");
	hbt.randy->reset(-time(NULL));
	hbt.ReadPR();
	printf("check in\n");
	hbt.CalcCF_MC();
	printf("check out\n");
	//hbt.CFArray[0][0][10]->Print();
	printf("check before averaging\n");
	hbt.AverageCF();
	return 0;
}
