#ifndef __HBT_BES_H__
#define __HBT_BES_H__
#include "commonutils.h"
#include "coral.h"
#include "parameterMap.h"
#include "randy.h"

using namespace std;

class CHBT_Part{
public:
	CHBT_Part();
	int ID;
	double mass,pt;
	vector<double> p;
	vector<double> xbj;
	void Print();
};

class CHBT_BES{
public:
	CparameterMap *parmap;
	CWaveFunction *wf;
	CRandy *randy;
	
	string INPUT_OSCAR_FILENAME;
	string OUTPUT_CF_FILENAME;
	double TAU_COMPARE,QINVTEST;
	
	// CF objects
	int Nqinv,IDA,IDB,NEVENTS_MAX,NMC;
	double DELqinv,qinvMAX,RANSEED;
	vector<double> CFqinv;
	
	// Part List
	vector<CHBT_Part *> parta;
	vector<CHBT_Part *> partb;
	
	CHBT_BES(string parsfilename);
	void ReadPR();
	double Getqinv(vector<double> &pa,vector<double> &pb);
	void CalcXR(CHBT_Part *partaa,CHBT_Part *partbb,vector<double> &x,double &r);
	void AddPart(vector<CHBT_Part *> &part,int &IDread,vector<double> &pread,vector<double> &xread);
	void CalcXBjPt(CHBT_Part *partaa,vector<double> &xread);
	void CalcCF();
	void CalcCF_MC();
};

#endif
