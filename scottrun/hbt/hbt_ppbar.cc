#include "hbt_bes.h"

using namespace std;

int main(){
	CHBT_BES hbt("parameters.txt");
	printf("howdy\n");
	hbt.ReadPR();
	hbt.CalcCF();
	return 0;
}
