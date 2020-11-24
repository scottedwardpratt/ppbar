#ifndef __HBT_BES_H__
#define __HBT_BES_H__
#include "commonutils.h"
#include "coral.h"
#include "parametermap.h"
#include "randy.h"

using namespace std;
class CF;
class CHBT_Part;
typedef multimap<double,CHBT_Part *> CHBT_PartMap;

class CHBT_Part{
public:
	CHBT_Part();
	double rdummy;
	int ID;
	double mass,pt,uperp;
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
	bool GAUSS;  // if true don't init stuff for reading in OSCAR output
	
	string INPUT_OSCAR_FILENAME;
	string INPUT_OSCAR_BASE_DIRECTORY;
	string RESULTS_DIR;
	int NRAP,NUPERP,NPHI;
	double MASSA,MASSB;
	// YMAX should be BIGGER Than DELRAP*NRAP (so that for y>YMAX it can never average with another particle to be inside rapidity window)
	double DELRAP,DELUPERP,DELPT,DELPHI,YMAX; 
	int INPUT_OSCAR_NRUNS;
	string OUTPUT_CF_FILENAME;
	string GITHOME_MSU;
	double TAU_COMPARE,UPERPTEST;
	
	// CF objects
	int IDA,IDB,QAB,NEVENTS_MAX,NEVENTS;
	long long NMC,NTRY;
	double RANSEED;

	// Part List
	//vector<CHBT_Part *> parta;
	//vector<CHBT_Part *> partb;
	vector<vector <CHBT_PartMap>> partmap;
	vector<vector<vector<int>>> Na,Nb;
	
	CHBT_BES(string parsfilename);
	void ReadPR();
	void ReadCFs();
	double Getqinv(vector<double> &pa,vector<double> &pb);
	void AddPart(int &IDread,vector<double> &pread,vector<double> &xread,double mass);
	void CalcXBjPt(CHBT_Part *partaa,vector<double> &xread);
	void CalcCF();
	void CalcCF_Gauss(double Rout,double Rside, double Rlong);
	void WriteCFs();
	CF *GetCF(CHBT_Part *parta,CHBT_Part *partb);
	void AverageCF();
	void CalcCoalescenceSpectra();
	void GetIrapIphiIuperp(double rap,double phi,double uperp,int &irap,int &iphi,int &iuperp);
	vector<vector<vector<CF *>>> CFArray;
	void WriteThetaPhiDists();
	CF *cfbar;
	CF *cfgauss;
};

class CF{
public:
	static int NQ,Nxyz,NSAMPLE_THETAPHI;
	static double DELQ;
	static double OUTSIDELONG_DIRECTION_CUT;  // cos(theta) must be > this value, theta is angle of q rel to axis
	static double OUTSIDELONG_Q_CUT;
	static double Dxyz;
	static bool USE_OUTSIDELONG_DIRECTION_CUT;
	static bool USE_OUTSIDELONG_Q_CUT;
	static double Rcoalescence;
	static CRandy *randy;
	static CHBT_BES *hbt;
	static bool COAL_USE_WF;
	static double COAL_DELR;
	static vector<double> psi_coal;
	void CalcCoalWF();
	double D3q;  // Used for weighting and coalescence
	long long int nincrement;
	double ncoalescence;
	double CoalescenceWeight(double r);
	CWaveFunction *wf;
	vector<double> cf_qinv,cf_qout,cf_qside,cf_qlong;
	vector<double> source_out,source_side,source_long;
	vector<double> norm_qinv,norm_qout,norm_qside,norm_qlong;
	vector<vector<double>> ThetaPhiDist;
	CF();
	void CalcXR(CHBT_Part *partaa,CHBT_Part *partbb,vector<double> &x,double &r);
	void Reset();
	void Normalize();
	void Print();
	void WriteCFs(string filename);
	void ReadCF(string filename);
	void WriteThetaPhiDist(string filename);
	void PrintSourceProjections();
	void Increment(CHBT_Part *parta,CHBT_Part *partb);
	void Increment(vector<double> &x);
};

#endif
