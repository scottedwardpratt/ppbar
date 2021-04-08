#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <boost/math/special_functions.hpp>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;
using namespace boost::math;

int main(int argc,char *argv[]){
  complex<double> qa,B;
	double qaR,qaI,qmaga,phi,phismallest=0.0,dphi=0.0001*PI;
	complex<double> ci(0.0,1.0);
	double smallestB2,B2;
	double ka,sigma;
	
	
	printf("Enter phi in degrees: ");
	scanf("%lf",&phi);
	phi=PI*phi/180.0;

	/*
	for(ka=0.1;ka<10;ka+=0.2){
		printf("--------- ka=%g __________\n",ka);
		for(qmaga=0.25*ka;qmaga<100;qmaga+=0.2*ka){
			qaR=qmaga*cos(phi);
			qaI=qmaga*sin(phi);
			qa=qaR+ci*qaI;
			B=qa*cos(qa)+ci*ka*sin(qa);
			B=B/(qa*cos(qa)-ci*ka*sin(qa));
			B*=exp(-2.0*ci*ka);
			B2=real(B*conj(B));
			printf("qmag=%6.3f: qa=(%g,%g),  B=(%g,%g), B^2=%g\n",qmaga,real(B),imag(B),real(qa),imag(qa),B2);
		}
	}
	*/
	
	printf("Enter qmaga: ");
	scanf("%lf",&qmaga);
	qaR=qmaga*cos(phi);
	qaI=qmaga*sin(phi);
	qa=qaR+ci*qaI;
	for(ka=0.2;ka<10;ka+=0.2){
		B=qa*cos(qa)+ci*ka*sin(qa);
		B=B/(qa*cos(qa)-ci*ka*sin(qa));
		B*=exp(-2.0*ci*ka);
		B2=real(B*conj(B));
		sigma=2.0*real(1.0-B)-real((1.0-B)*conj(1.0-B));
		sigma=PI*sigma/(ka*ka);
		printf("ka=%7.3f:  B=(%10.5f,%10.5f), B^2=%10.5f, sigma/a^2=%10.5f, (sigma/a^2)*pow(ka,0.7)=%10.5f\n",ka,real(B),imag(B),B2,sigma,sigma*ka);
	}
	
	
	
	/*
	for(ka=0.25;ka<10;ka*=2){
		printf("------------ ka=%g ---------------\n",ka);
		for(qmaga=0.2*ka;qmaga<20*ka;qmaga+=0.2*ka){
			smallestB2=1000000000.0;
			for(phi=-PI+dphi;phi<PI;phi+=dphi){
				qaR=qmaga*cos(phi);
				qaI=qmaga*sin(phi);
				qa=qaR+ci*qaI;
				B=qa*cos(qa)+ci*ka*sin(qa);
				B=B/(qa*cos(qa)-ci*ka*sin(qa));
				B*=exp(-2.0*ci*ka);
				B2=real(B*conj(B));
				if(B2<smallestB2){
					smallestB2=B2;
					phismallest=phi;
				}
			}
			printf("qmaga=%10.6f: smallestB2=%10.5f, phi=%8.4f\n",qmaga,smallestB2,phismallest*180.0/PI);
		}
	}
	*/
	
  return 0;
}





