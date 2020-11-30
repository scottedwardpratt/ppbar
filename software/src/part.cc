#ifndef __HBT_BES_PART_CC__
#define __HBT_BES_PART_CC__
#include <cstring>
#include "commondefs.h"
#include "coral.h"
#include "hbt_bes.h"

CHBT_Part::CHBT_Part(){
	p.resize(4);
	x.resize(4);
}

void CHBT_Part::Print(){
	printf("ID=%d, mass=%g\n",ID,mass);
	printf("x=(%g,%g,%g,%g), p=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3],p[0],p[1],p[2],p[3]);
}

#endif