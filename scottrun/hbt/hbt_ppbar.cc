#include "hbt_bes.h"

using namespace std;

int main(){
	CHBT_BES hbt("parameters.txt");
	hbt.randy->reset(-time(NULL));
	hbt.ReadPR();
	hbt.CalcCF_MC();
	return 0;
}
