#include "hbt_bes.h"

using namespace std;

int main(){
	CHBT_BES hbt("parameters.txt");
	printf("howdy\n");
	hbt.randy->reset(-time(NULL));
	hbt.ReadPR();
	hbt.CalcCF_MC();
	return 0;
}
