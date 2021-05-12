#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <boost/math/special_functions.hpp>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;
const double MPROTON=938.28;

using namespace std;
using namespace boost::math;
double GetK(double plab){
	double mass=939.0,e1,s,k;
	e1=sqrt(mass*mass+plab*plab);
	s=(e1+mass)*(e1+mass)-plab*plab;
	e1=0.5*sqrt(s);
	k=sqrt(e1*e1-mass*mass);
	return k;
}

void GetBessel(unsigned int ellmax,complex<double> x,vector<complex<double>> &jl,vector<complex<double>> &jlprime){
	unsigned int ell;
	complex<double> xx,y0,y1,ci(0.0,1.0);
	
	if(jl.size()<ellmax+1){
		jl.resize(ellmax+1);
		jlprime.resize(ellmax+1);
	}
	jl[0]=sin(x)/x;
	jl[1]=sin(x)/(x*x)-cos(x)/x;
	jlprime[0]=-jl[1];
	jlprime[1]=jl[0]-2.0*jl[1]/x;
	for(ell=2;ell<=ellmax;ell++){
		jl[ell]=-jl[ell-2]+(2.0*ell-1.0)*jl[ell-1]/x;
		jlprime[ell]=jl[ell-1]-(ell+1.0)*jl[ell]/x;
	}
}

void GetHankel(unsigned int ellmax,double x,vector<complex<double>> &hl,vector<complex<double>> &hlprime){
	unsigned int ell;
	complex<double> ci(0.0,1.0);
	double j0,j1,y0,y1;
	
	if(hl.size()<ellmax+1){
		hl.resize(ellmax+1);
		hlprime.resize(ellmax+1);
	}
	j0=sin(x)/x;
	j1=sin(x)/(x*x)-cos(x)/x;
	y0=-cos(x)/x;
	y1=-cos(x)/(x*x)-sin(x)/x;
	hl[0]=j0+ci*y0;
	hl[1]=j1+ci*y1;
	hlprime[0]=-hl[1];
	hlprime[1]=hl[0]-2.0*hl[1]/x;
	for(ell=2;ell<=ellmax;ell++){
		hl[ell]=-hl[ell-2]+(2.0*ell-1.0)*hl[ell-1]/x;
		hlprime[ell]=hl[ell-1]-(ell+1.0)*hl[ell]/x;
	}
}

double GetSigma(double VR,double VI,double a,double plab){
	double qmaga,qmag2,ka,delsigma,sigma;
	complex<double> qa,B,numer,denom,ci(0.0,1.0);
	vector<complex<double>> jl,hl,jlprime,hlprime;
	double mu=MPROTON*0.5,k,phi;
	unsigned int ell,ellmax=10;
	k=GetK(plab);
	ka=k*a/HBARC;
	qmag2=2.0*mu*VI+(k*k-2.0*mu*VR);
	qmaga=sqrt(qmag2)*a/HBARC;
	phi=0.5*atan2(2.0*mu*VI,k*k-2.0*mu*VR);
	if(sin(phi)<0.0){
		printf("what???? sin(phi)=%g, phi=%g\n",sin(phi),phi*180.0/PI);
		exit(1);
	}
	
	qa=qmaga*cos(phi)+ci*qmaga*sin(phi);
	sigma=0.0;
	GetBessel(ellmax,qa,jl,jlprime);
	GetHankel(ellmax,ka,hl,hlprime);
	for(ell=0;ell<=ellmax;ell++){
		denom=hl[ell]*jlprime[ell]*qa-ka*jl[ell]*hlprime[ell];
		numer=conj(hlprime[ell]*ka)*jl[ell]-conj(hl[ell])*jlprime[ell]*qa;
		//printf("ell=%d, ka=%g, qmaga=%g, qa=(%g,%g), denom=(%g,%g)\n",ell,ka,qmaga,real(qa),imag(qa),real(denom),imag(denom));
		B=numer/denom;
		//B=B*pow(ci,ell);
		delsigma=2.0*real(1.0-B)-real((1.0-B)*conj(1.0-B));
		if(delsigma<-1.0E-3){
			printf("delsigma <0!!!\n");
			printf("ell=%d, plab=%g, delsigma=%g\n",ell,plab,delsigma);
			printf("B=(%g,%g)\n",real(B),imag(B));
			printf("ka=%g, qa=(%g,%g), phi=%g\n",ka,real(qa),imag(qa),phi*180.0/PI);
			printf("VR=%g, VI=%g\n",VR,VI);
			exit(1);
		}
		sigma+=(2.0*ell+1)*delsigma;
	}
	sigma=a*a*PI*sigma/(ka*ka);
	return sigma;
}

double GetChiSquare(double VR,double VI,double a){
	double plab,sigma,chisquare=0.0,pmax=1000.0,formula;
	for(plab=150;plab<pmax;plab+=50){
		sigma=GetSigma(VR,VI,a,plab);
		formula=6.7*pow(plab/1000.0,-0.7);
		chisquare+=pow(formula-sigma,2.0)/formula;
	}
	return chisquare;
}

int main(int argc,char *argv[]){
	double chisquare,bestchisquare=1.0E10,VR,VI,a,besta=1.0E10,bestVR=1.0E10,bestVI=1000000.0,formula,plab,sigma;
	for(a=1.81;a<1.85;a+=0.001){
		for(VR=-9.5;VR<-8.0;VR+=0.01){
			for(VI=24.5;VI<25.5;VI+=0.01){
				chisquare=GetChiSquare(VR,VI,a);
				if(chisquare<bestchisquare){
					besta=a;
					bestVR=VR;
					bestVI=VI;
					bestchisquare=chisquare;
					printf("VR=%g, VI=%g, a=%g, chisquare=%g\n",VR,VI,a,chisquare);
				}
			}
		}
	}
	printf("best a=%g fm, best VR=%g, best VI=%g, best chi^2=%g\n",besta,bestVR,bestVI,bestchisquare);
	
	a=besta;
	VR=bestVR;
	VI=bestVI;
	for(plab=50;plab<1500;plab+=50){
		formula=6.7*pow(plab/1000.0,-0.7);
		sigma=GetSigma(VR,VI,a,plab);
		printf("%6.1f %6.1f %8.4f %8.4f\n",plab,GetK(plab),sigma,formula);
	}
	
	return 0;
}





