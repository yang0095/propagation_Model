// multDlg.h : header file
//

#if !defined(AFX_MULTDLG_H__C5EB711A_71A4_4692_9558_7D72CC108275__INCLUDED_)
#define AFX_MULTDLG_H__C5EB711A_71A4_4692_9558_7D72CC108275__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CMultDlg dialog

typedef struct {double real, imag;} dxcomplex;


class CMultDlg : public CDialog
{
// Construction
public:
	int GetFileData(CString sPath, double** fData1, double** fData2);
	void CalAttIo(double f, double theta, double *attio);
	void ocdiv(dxcomplex *a, dxcomplex *b);
	void ontrt(dxcomplex *a, int n);
	void CalAttMult(double theta, int season, double dsa, double hs, double ha, int ipol, double pout, double f, double* dt, double* damp, double* AMout,double* Lr);
	void newCalAttMS(double rou0, double Ts, int season, double hs, double lat, double theta, double f, double p0, double d, double yita, double *attsint);
	void CalAttMS(int season, double Ts, double rou0, double lat, double theta, double f, double pout, double d, double yita, double *AMSout);
	void SintAtt(int season, double Ts, double rou0, double lat, double theta, double f, double p0, double d, double yita, double *attsint, double *newp);
	void CalAttGas(double ele,double h_s,double h_areo, double f, double T, double P, double rou0, double *attgas);
	void spec_gasatt(double freq, double T, double Rou, double P, double *gama);
	void CalAttCloud(double Ts, double Ps, double rou0, double rain_h,double plat_h1,double plat_h2,double f,double theta0,double LW,double *attcloud);
	void CalAttRain(double lats, double theta0, double ipol, double f, double hs, double hr,double ha, double R001, double p, double *attrain);
	void CalFrady(double N[], double h[], double ht, double f, double theta, double fai, double *Frady, double *attfrady, double *XPD);
	void CalSeSan(float N[], float h[], float ht, float f, float df, float theta, float *ShiY, float *XiangW, float *Pe);
    float TEC(float N[], float h[], float ht, float theta);
	void CalShiY(float N[], float h[], float ht, float f, float theta, float *ShiY);
    void CalIoSin(double l, double fai, double fs, double hs, double f, double theta, double *Ais);
	CMultDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CMultDlg)
	enum { IDD = IDD_MULT_DIALOG };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMultDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CMultDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnButton1();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MULTDLG_H__C5EB711A_71A4_4692_9558_7D72CC108275__INCLUDED_)
