#ifndef __HBT_BES_H__
#define __HBT_BES_H__
#include "commonutils.h"
#include "coral.h"
#include "parametermap.h"
#include "randy.h"

using namespace std;
class CF;

class CHBT_Part{
public:
	CHBT_Part();
	double rdummy;
	int ID;
	double mass,pt;
	vector<double> p;
	vector<double> x;
	double rap0,phi0; // before particles boosted to rap=0 and phi=0.
	void Print();
};

class CHBT_BES{
public:
	CparameterMap *parmap;
	CWaveFunction *wf;
	CRandy *randy;
	
	string INPUT_OSCAR_FILENAME;
	string INPUT_OSCAR_BASE_DIRECTORY;
	string RESULTS_DIR;
	int NRAP,NPT,NPHI;
	// YMAX should be BIGGER Than DELRAP*NRAP (so that for y>YMAX it can never average with another particle to be inside rapidity window)
	double DELRAP,DELPT,YMAX;
	int INPUT_OSCAR_NRUNS;
	string OUTPUT_CF_FILENAME;
	string GITHOME_MSU;
	double TAU_COMPARE,QINVTEST;
	
	// CF objects
	int IDA,IDB,NEVENTS_MAX,NMC;
	double RANSEED;
	
	// Part List
	vector<CHBT_Part *> parta;
	vector<CHBT_Part *> partb;
	
	CHBT_BES(string parsfilename);
	void ReadPR();
	double Getqinv(vector<double> &pa,vector<double> &pb);
	void AddPart(vector<CHBT_Part *> &part,int &IDread,vector<double> &pread,vector<double> &xread);
	void CalcXBjPt(CHBT_Part *partaa,vector<double> &xread);
	void CalcCF();
	void CalcCF_MC();
	void WriteCFs();
	CF *GetCF(CHBT_Part *parta,CHBT_Part *partb);
	void AverageCF();
	vector<vector<vector<CF *>>> CFArray;
	void WriteThetaPhiDists();
	CF *cfbar;
};

class CF{
public:
	static int NQ,Nxyz;
	static double DELQ;
	static double OUTSIDELONG_DIRECTION_CUT;  // cos(theta) must be > this value, theta is angle of q rel to axis
	static double OUTSIDELONG_Q_CUT;
	static double Dxyz;
	static bool USE_OUTSIDELONG_DIRECTION_CUT;
	static bool USE_OUTSIDELONG_Q_CUT;
	static CRandy *randy;
	static CHBT_BES *hbt;
	long long int nincrement;
	CWaveFunction *wf;
	vector<double> cf_qinv,cf_qout,cf_qside,cf_qlong;
	vector<double> source_out,source_side,source_long;
	vector<long long int> norm_qinv,norm_qout,norm_qside,norm_qlong;
	vector<vector<double>> ThetaPhiDist;
	CF();
	void CalcXR(CHBT_Part *partaa,CHBT_Part *partbb,vector<double> &x,double &r);
	void Reset();
	void Normalize();
	void Print();
	void WriteCFs(string filename);
	void WriteThetaPhiDist(string filename);
	void PrintSourceProjections();
	void Increment(CHBT_Part *parta,CHBT_Part *partb);
};

#endif
