// P1812_Fun.h: interface for the CP1812_Fun class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_P1812_FUN_H__72A1B6E6_62E7_462E_9537_E3AF2AE56911__INCLUDED_)
#define AFX_P1812_FUN_H__72A1B6E6_62E7_462E_9537_E3AF2AE56911__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CP1812_Fun  
{
public:
	CP1812_Fun();
	virtual ~CP1812_Fun();

	bool P1812_6_ITUwb(double freq, double p, double pL, double latMid, double htg, double hrg, int pol, double erp, int n, double *dn, double *hn,double *Rn, int *climn, double dN, double N0,
		double &Lb,double &Ep, CString &errorMes); 
	bool P1812_6_ITUwb(double freq, double p, double pL, double latMid, double htg, double hrg, int pol, double erp, int n, double *dn, double *hn,double *Rn, int *climn, double dN, double N0,
							   double &ae, bool &bpath, double &thetat, double &thetar, double &theta, double &hm, double &Lb,double &Ep, CString &errorMes); 	
	void RadioClimZone1(int n,double *dn,double *hn,int *ngc,int *climn);

	//TEMP
	double I(double x);


private:
	void ParaInit(double latMid,int n,double *dn, int *climn, double dN, double &ae,double &b0,double &tau,double &w);
	void GetDtmDlm(int n,double *dn,int *climn,double &dtm, double &dlm);
	double GetW(int n, double *dn, int *climn);
	double GetAe(double dN);
	void GetB0(double dlm, double dtm, double latMid,double &b0,double &tau);
	void CompLossTime(double freq,double p, int pol,double htg,double hrg,int n,double *dn, double *hn,double *gn, int *climn,double N0,double ae,double b0,double tau,double w,
		double &lb0p,double &lbc);
	double CompLossLoc(double freq, double pL, double h, double R, double lb0p,double lbc);
	void PathAnalyze(double freq,double hts,double hrs,double ae,int n,double *dn, double *hn,bool &bpath, int &ilt, int &ilr, double &thetat, double &thetar, double &theta);
	int Im50(double freq,double h1,double h2, int n, double *dn, double *hn, double ae);
	double LosFreespace(double freq,double d,double p,double dlt,double dlr);
	void Ld1812_3(double freq,double htc,double hrc,int pol,double p,double b0,double ae,double w,int n,double *dn, double *hn, double *gn,double &ld50, double &ldp);
	double LdDeltaBullington(double freq,double htc,double hrc,int pol,double ap,double w,int n,double *dn,double *hn,double *gn);
	double LdBullington(double freq,double htc,double hrc,int n,double *dn,double *gn,double ae);
	double LdSphericalEarth(double freq,int pol,double htesph,double hresph, double d,double ae,double w);
	double GetFi(double p, double b0);
	void DeriSmoothEarthsurface(int n, double *dn, double *hn, double &hst, double &hsr);
	void DeriDuctPar(int n, double *dn,double *hn, double hts, double hrs, int ilt, int ilr, double hst, double hsr,double &hm, double &hte,double &hre);
	void DeriDiffPar(double hst,double hsr,double htc,double hrc,int n,double *dn,double *hn,double &hstd, double &hsrd);
	double LosDuct(double freq,double p,double b0,double hts,double hrs,int n, double *dn,double *hn, int *climn,double ae,double w,double tau, int ilt, int ilr, double thetat, double thetar);
	void GetDctDcr(int n, double *dn, int *climn, double &dct, double &dcr);
	double Stim(double htc,int n, double *dn, double *gn, double ae);
	double Srim(double hrc, int n, double *dn, double *gn, double ae);
	double Vmax(double freq,double h1,double h2, int n, double *dn, double *gn, double ae);
	double LosKnifeEdge(double v);
	double LdFirstTermComp(double freq, double d,double htesph,double hresph, double adft,double w,int pol);
	double LdFirstTerm(double freq, double d,double htesph,double hresph,double eps,double sgm, double adft,int pol);
	double LosTropoScatter(double freq, double p, double d, double theta, double N0);
	double GetFjk(double Theta, double xi, double theta);
	//double I(double x);
	double L2E(double freq, double lb);
	double T(double x);
	double Xi(double x);
	
};

#endif // !defined(AFX_P1812_FUN_H__72A1B6E6_62E7_462E_9537_E3AF2AE56911__INCLUDED_)
