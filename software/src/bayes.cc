#ifndef __HBT_BES_BAYES_CC__
#define __HBT_BES_BAYES_CC__
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

double CHBT_BES::GetChiSquare(CF *cfa,CF *cfb){
	double chisquare=0.0;
	int iq;
	double diff;
	for(iq=0;iq<CF::NQ;iq++){
		diff=cfa->cf_qout[iq]-cfb->cf_qout[iq];
		chisquare+=diff*diff;
		diff=cfa->cf_qside[iq]-cfb->cf_qside[iq];
		chisquare+=diff*diff;
		diff=cfa->cf_qlong[iq]-cfb->cf_qlong[iq];
		chisquare+=diff*diff;
	}
	return chisquare;
}


#endif