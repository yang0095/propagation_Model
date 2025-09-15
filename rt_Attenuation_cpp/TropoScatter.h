// TropoScatter.h: interface for the CTropoScatter class.
//
//////////////////////////////////////////////////////////////////////


#if !defined(AFX_TROPOSCATTER_H__6BC03E48_E74A_4D8F_8DF2_832E7E4D6511__INCLUDED_)
#define AFX_TROPOSCATTER_H__6BC03E48_E74A_4D8F_8DF2_832E7E4D6511__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000




class CTropoScatter  
{
public:
	void TropoScatCal(double freq, double hts, double hrs, double d, double p,  double Gt, double Gr,double dN, int climate, double &los, bool &ipath);
	CTropoScatter();
	virtual ~CTropoScatter();

private:	
	double Xi(double x);
	double T(double x);
	double Y90(double freq, double h, int climate, double ae, double Theta0);
	double Qi(double x);
	void GetPara(int climate, double &M, double &gamma);
	void LosTroposcatter(double freq, double d, double p,  double Theta0,double Gt, double Gr,double ae,
		int climate,double &los, double &sigma);
	double GetAe(double dN);
	void Interp1(double *x, double *y, int n, double *x1, double *y1, int n1, CString method);
	double Linear(double x1, double y1, double x2, double y2, double x0);


};

#endif // !defined(AFX_TROPOSCATTER_H__6BC03E48_E74A_4D8F_8DF2_832E7E4D6511__INCLUDED_)
