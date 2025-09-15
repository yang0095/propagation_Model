// multDlg.cpp : implementation file
//

#include "stdafx.h"
#include "mult.h"
#include "multDlg.h"
#include "math.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMultDlg dialog

CMultDlg::CMultDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CMultDlg::IDD, pParent)
{

	//{{AFX_DATA_INIT(CMultDlg)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMultDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CMultDlg)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CMultDlg, CDialog)
	//{{AFX_MSG_MAP(CMultDlg)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON1, OnButton1)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMultDlg message handlers

BOOL CMultDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	
	// TODO: Add extra initialization here
	
	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CMultDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CMultDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CMultDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CMultDlg::OnButton1() 

{

	CMultDlg Dlg;
	double re=6371000;//地球半径，m
    double pi=3.1415926;


///////////////////////////////////////用户输入参数//////////////////////////////////////    
	double lats=36.1;//地面节点纬度，（0，60）度
	double lons=120.3;//地面节点经度，（0，180）度
	double hs=7;//地面节点离地高度（0，3000）m


	double lata=36.1;//空中节点纬度，（0，60）度25.8 25.68 25.1 24.95 24.72 24.4 24.2 24.04 23.93 23.905 23.8
	double lona=120.3;//空中节点经度，（0，180）度
	double ha=40000;//空中节点离地高度（3000,10000）m

	int season=1;//1冬季，2春季，3夏季，4秋季

	double pout=99;//衰减未被超过的时间概率[80，99]
	double f=10;//频率, (15,30)GHz 0.45 1.3
	double d=1;//天线直径, (0.1,5)m
	double yita=0.5;//天线效率 ,  (0.1,1)
	int ipol=0;//天线极化方式，0：水平极化，1：垂直极化，2：圆极化

///////////////////////需要从数据库中读取的中间输入参数/////////////////////////////////////////////

	double Ts=-2.1+273.15;//地表面温度（K），基于季节season和地表面经纬度lons，lats从数据库中获取
	double Ps=1014;//地表面大气压强（hpa），基于季节season和地表面经纬度lons，lats从数据库中获取
	double rou0=1.9;//水汽密度(g/m^3),基于概率pout和地表面经纬度lons，lats从数据库中获取

	double water=3;//液态水含量(kg/m^2)基于概率pout和地表面经纬度lons，lats从数据库中获取99%:3  50%:0
	double hr=5000;//雨顶高度hr（m）,基于地表面经纬度lons，lats从数据库中获取
	double R001=100;//分钟降雨率(mm/h),基于地表面经纬度lons，lats从数据库中获取



////////////////计算地心角fai(rad)，路径距离dsa(km),真实仰角theta0(度),视在仰角theta(度),////////////////////////////////////////////////////////////////

	bool iout=false;//目标在视距之外标志，若为true，判为出错，弹出提示框“目标在视距之外！”
	double attmulsin=0;
    double dlat=0;
	double dlon=0;
	double a1=0;
	double fai=0;
	double dsa=0;
	double theta=0;
	double r1=0;
	double r2=0;
	double theta0=0;
	double tao=0;
	double thetam=0;
	double taos=0;
	double dt=0;
	double damp=0;
	double attmul=0;
	double attcloud=0;
	double attio=0;
	double gama=0;
	double attgas=0;
	double attrain=0;
	double attiscin=0;

	double attall=0;
	double Lr=0;

	double Frady=0;
	double attfrady=0;
	double XPD=0;



	CString D = "F:\\ion.txt";//电子浓度剖面地址

	int nD = 0;        //电子浓度剖面层数
	double* D1 = NULL; //电子浓度剖面第1列数组（高度km）
	double* D2 = NULL; //电子浓度剖面第2列数组（电子浓度1/cm3）
	double belev = 0; 
	

	nD = GetFileData(D, &D1, &D2);


    r1=re+hs;
	r2=re+ha;
	dlat=lats-lata;
	dlon=lons-lona;
	a1=pow(sin(pi/180.0*dlat/2),2)+cos(lata*pi/180.0)*cos(lats*pi/180.0)*pow(sin(pi/180.0*dlon/2),2);

//	fai=2*atan(sqrt(a1)/sqrt(1-a1));//地心角，rad
//	dsa=sqrt(r1*r1+r2*r2-2*r1*r2*cos(fai))/1000;	//路径距离,km

	dsa=200;
	fai=acos((r1*r1+r2*r2-dsa*dsa*1000000)/(2*r1*r2));
	
	//(计算内容三时，地心角由输入路径距离dsa计算：fai=acos((r1*r1+r2*r2-dsa*dsa)/(2*r1*r2));

	
////////////////////////////////////////////////////
	int nd=33;
	double* Adsa = NULL;
	Adsa = new double[nd];
	memset(Adsa, 0, sizeof(double) * nd); 

	double* Aattall = NULL;
	Aattall = new double [nd];


	double* Atheta0 = NULL;
	Atheta0 = new double [nd];

	double* Atheta = NULL;
	Atheta = new double [nd];

	double* attMS = NULL;
	attMS = new double [nd];

	double* attCL = NULL;
	attCL = new double [nd];

	double* attGAS = NULL;
	attGAS = new double [nd];


	double* attRN = NULL;
	attRN = new double [nd];

	double* attIO = NULL;
	attIO = new double [nd];


	double* attFRADY = NULL;
	attFRADY = new double [nd];

	double* attIS = NULL;
	attIS = new double [nd];


	double* FRADY = NULL;
	FRADY = new double [nd];

	Atheta[0]=0.5;//距离数列赋值1~450km，步进1km
	for(int ii=1; ii<nd;ii++)
	{
		if(Atheta[ii-1]<3)
			Atheta[ii]=Atheta[ii-1]+0.5;
		else
			Atheta[ii]=Atheta[ii-1]+2;
	}


	CString str;
	CStdioFile fileout;	


	CString sPath1 = "F:\\attALL.txt";//
	fileout.Open(sPath1,CStdioFile::modeReadWrite);	


double ffai=0;//方位角
for( ii=0; ii<nd;ii++)
{
	theta=Atheta[ii];
	if(theta==90)
	{
		theta=89.9;
	}
	
	Dlg.CalIoSin(lons, lats, lons, hs/1000, f, theta, &attiscin);		
	Dlg.CalFrady(D2, D1, ha, f*pow(10,9), theta, ffai, &Frady, &attfrady, &XPD);
	Dlg.newCalAttMS(rou0, Ts, season, hs, lats, theta, f, pout, d, yita, &attmulsin);
	Dlg.CalAttGas(theta, hs, ha, f, Ts, Ps, rou0, &attgas);
	Dlg.CalAttRain(lats, theta, ipol, f, hs, hr,ha, R001, pout, &attrain);
	Dlg.CalAttCloud(Ts, Ps, rou0, hr,hs,ha,f,theta,water,&attcloud);///
//	Dlg.CalAttMult(theta, season, dsa/1000, hs, ha, ipol, pout, f, &dt, &damp, &attmul,&Lr);	
//	Dlg.CalAttIo(f,theta,&attio);
//	double Lbf=20*log10(f*1000)+20*log10(dsa)-27.56;//自由空间衰减，dB
//   Aattall[ii]=attfrady+attio;
//   if(pout>=50)

//	else
//	{
//		if(pow(attcloud+attrain,2)-pow(attmulsin,2)>0)
//			Aattall[ii]=attgas+sqrt(fabs(pow(attcloud+attrain,2)-pow(attmulsin,2)));
//		else
//			Aattall[ii]=attgas-sqrt(fabs(pow(attcloud+attrain,2)-pow(attmulsin,2)));
//	}
		
//	Aattall[-1]=Aattall[0];
//	if(Aattall[ii]<Aattall[ii-1])
//			Aattall[ii]=Aattall[ii-1];

	attMS[ii]=attmulsin;
	attGAS[ii]=attgas;
	attRN[ii]=attrain;
	attCL[ii]=attcloud;
//	attIO[ii]=attio;
	attFRADY[ii]=attfrady;
	attIS[ii]=attiscin;

//	FRADY[ii]=Frady;

//	if(ii>1)
//	{
//		if(attMS[ii]>attMS[ii-1])
//		{
//			attMS[ii]=attMS[ii-1];
//		}
//		if(attRN[ii]>attRN[ii-1])
//		{
//			attRN[ii]=attRN[ii-1];
//		}
//	}	

		if(hs<3000) //地空链路
			Aattall[ii]=attGAS[ii]+sqrt(pow(attCL[ii]+attRN[ii],2)+pow(attMS[ii],2));//+attIS[ii]+attFRADY[ii]
		else//空空链路
			Aattall[ii]=attGAS[ii]+attCL[ii];
//	Atheta0[ii]=theta0;
//	Atheta[ii]=theta;
//str.Format("%f %f\r\n",Atheta[ii],attGAS[ii]);// , Atheta0[ii],
	str.Format("%f %f %f %f %f %f \r\n",Atheta[ii], attGAS[ii],attCL[ii],attRN[ii],attMS[ii],Aattall[ii]);//, Atheta0[ii],
//		str.Format("%f %f \r\n",Atheta[ii], FRADY[ii]);// , Atheta0[ii],

	fileout.WriteString(str);
}
	




////////////////////////////计算多径与闪烁衰落深度attmulsin(dB)////////////////////////////////////////////////////////////////
//	Dlg.CalAttMS(season, lats, theta, f, pout, d, yita, &attmulsin);

	
//////计算多径小尺度衰落深度attmul(dB)，内容四多径信道模型的相对时延dt(ns)，相对幅度damp(dB),多径损耗Lr(dB)////////////////////////////////////////////////////////////////

//	Dlg.CalAttMult(season, dsa, hs, ha, ipol, pout, f, &dt, &damp, &attmul,&Lr);	


////////////////////////////计算云雾衰减attcloud(dB)////////////////////////////////////////////////////////////////


//	Dlg.CalAttCloud(Ts,f,theta,water,&attcloud);

////////////////////////////计算大气吸收衰减attgas(dB)////////////////////////////////////////////////////////////////


//	Dlg.CalAttGas(theta, hs/1000, ha/1000, f, Ts, Ps, rou0, &attgas);

	
	
////////////////////////////计算降雨衰减attrain(dB)////////////////////////////////////////////////////////////////
	
//	Dlg.CalAttRain(lats, theta, ipol, f, hs/1000, hr,ha/1000, R001, 100-pout, &attrain);


///////////////////////计算总衰减attall(dB)/////////////////////////////////
//	attall=attgas+sqrt(pow(attcloud+attrain,2)+pow(attmulsin,2));



//////////////////////////////////////////////////////////////////////////////////////////////
/*	int np=101;//时间概率
	double* th = NULL;
	th = new double[np];
	memset(th, 0, sizeof(double) * np); 

	double* attMS = NULL;
	attMS = new double [np];


	th[0]=0;
	for(ii=1; ii<np;ii++)
	{
		th[ii]=th[ii-1]+0.1;
	}
	


	CString str;
	CStdioFile fileout;	


	CString sPath1 = "G:\\attMS.txt";//
	fileout.Open(sPath1,CStdioFile::modeReadWrite);	


	for(int ip=0;ip<np;ip++)
	{
		Dlg.CalAttMS(season, lats,th[ip], f, pout, d, yita, &attMS[ip]);
		str.Format("%f %f\r\n", th[ip],attMS[ip]);

		fileout.WriteString(str);

	}	





	for(int i=1;i<np;i++)
	{
		if(attMS[i]>attMS[i-1])
			attMS[i]=attMS[i-1];
	}

*/





/*	CString sPath2 = "G:\\孙\\29suo201311\\40to99\\data\\attM40.txt";//输出数据文件：第一列p，第二列attMS
	fileout.Open(sPath2,CStdioFile::modeReadWrite);	


for(ifr=0; ifr<nf; ifr++)

{
	for(int ip=0;ip<np;ip++)
	{

		str.Format("%f\r\n", attM[ip][ifr]);

		fileout.WriteString(str);

	}
}

	CString sPath3 = "G:\\孙\\29suo201311\\40to99\\data\\attGAS40.txt";//输出数据文件：第一列p，第二列attMS
	fileout.Open(sPath3,CStdioFile::modeReadWrite);	


for(ifr=0; ifr<nf; ifr++)

{
	for(int ip=0;ip<np;ip++)
	{

		str.Format("%f\r\n", attGAS[ip][ifr]);

		fileout.WriteString(str);

	}
}


	CString sPath4 = "G:\\孙\\29suo201311\\40to99\\data\\attRN40.txt";//输出数据文件：第一列p，第二列attMS
	fileout.Open(sPath4,CStdioFile::modeReadWrite);	


for(ifr=0; ifr<nf; ifr++)

{
	for(int ip=0;ip<np;ip++)
	{
		str.Format("%f\r\n", attRN[ip][ifr]);

		fileout.WriteString(str);

	}
}

	CString sPath5 = "G:\\孙\\29suo201311\\40to99\\data\\attA40.txt";//输出数据文件：第一列p，第二列attMS
	fileout.Open(sPath5,CStdioFile::modeReadWrite);	


for(ifr=0; ifr<nf; ifr++)

{
	for(int ip=0;ip<np;ip++)
	{
		str.Format("%f\r\n", attA[ip][ifr]);

		fileout.WriteString(str);

	}
}

	fileout.Close();



	/*double* th = NULL;//时间概率
	th = new double[np];
	memset(th, 0, sizeof(double) * np); 

	th[0]=0;
	for(int ii=1; ii<np;ii++)
	{
		th[ii]=th[ii-1]+0.1;
	}*/
	

	//double thout=0;
	

//	for(i=1;i<np;i++)
//	{
//		if(attM[i]>attM[i-1])
//			attM[i]=attM[i-1];
//	}

}

///////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////

///////////////////////所有传播效应的计算函数：////////////////////////////////////////////////////////////////////


///////////////////////大气吸收计算函数CalAttGas////////////////////////////////////////////////////////////////////


void CMultDlg::CalAttGas(double ele,double h_s,double h_areo, double f, double T, double P, double rou0, double *attgas)
{
    //输入参数：真实仰角ele（度），离地高度h_s（km）、飞行器高度h_areo（km）、频率f（GHz）、温度T（k）、基于概率得到相应水汽密度rou0(g/m^3)、水汽柱积分含量
	//输出参数：大气衰减
	//ITU-R P.676-9 斜路径大气衰减的精确计算方法
	//定义变量
	CMultDlg Dlg;
	h_s=h_s/1000;
	h_areo=h_areo/1000;
	double gama=0;
	double pi=3.1415926;
	double r0=6371;
	double e;//定义水汽压
	int jj;//分层数目
	double delta1[1000];//用于确定分层数目
	ele=ele*pi/180;//仰角转化为弧度

	//依据飞行器高度确定分层数目
	double sum_delta=h_s;//定义总高度，不能超过飞行器高度
	for(int i=0;i<1000;i++)
	{
		delta1[i]=0.0001*exp(i*1E-2);
		for (int j=0;j<i+1;j++)
		{
			sum_delta=sum_delta+delta1[j];
		}
		if (sum_delta>=h_areo)
		{
			jj=i;
			break;
		}
		else
		{
			sum_delta=h_s;
		}
	}
	//定义中间变量
	
	double a[1000];
	double alfa[1000];
	double beta[1000];
	double r_abs;
	double delta[1000];
	double r[1000];
	double n[1000];
	double Th;
	double rou_h;
	double Ph;	
	double N;
	double TotalAtt=0;
	double temp_grad=-6.5;//11km以内温度线性递减
	//////	
	////////////
	//第一层初值
	r[0]=r0+h_s;//第一层地球中心到第一层的高度
	beta[0]=pi/2-ele;

	for (i=0;i<jj+1;i++)
	{
		delta[i]=0.0001*exp(i*1E-2);
		if(i==0)//分层第一层
		{		
			a[i] = -r[i]*cos(beta[i])+0.5*pow(4*pow(r[i],2.0)*pow(cos(beta[i]),2)+8.0*r[i]*delta[i]+4.0*pow(delta[i],2),0.5);
			alfa[i]=pi-acos(-(pow(a[i],2)+2*r[i]*delta[i]+pow(delta[i],2))/(2*a[i]*r[i]+2*a[i]*delta[i]));
			r_abs=r[i]-r0;//每层离地高度
			Th=T-6.5*r_abs;//离地高度为r_abs时的温度,压强和水汽密度
			Ph=P*pow(T/(T-6.5*r_abs),-34.163/6.5);
			rou_h=rou0*exp(-r_abs/2.0);
			e = rou_h*Th/216.7;
			N = 77.6/Th*(Ph+4810.0*e/Th);
			n[i] = 1+N*1E-6;			
			Dlg.spec_gasatt(f,Th, rou_h, Ph, &gama);
			TotalAtt = a[i]*gama+TotalAtt;
		}
		else
		{
			r[i] = r[i-1]+delta[i-1];
			r_abs = r[i]-r0;
			Th=T-6.5*r_abs;
			Ph=P*pow(T/(T-6.5*r_abs),-34.163/6.5);
			rou_h=rou0*exp(-r_abs*0.5);
			e = rou_h*Th/216.7;
			N = 77.6/Th*(Ph+4810.0*e/Th);
			n[i] = 1+N*1E-6;	
			beta[i] = asin(n[i-1]/n[i]*sin(alfa[i-1]));
			a[i]=-r[i]*cos(beta[i])+0.5*pow(4.0*pow(r[i],2)*pow(cos(beta[i]),2)+8.0*r[i]*delta[i]+4.0*pow(delta[i],2),0.5);
			alfa[i]=pi-acos(-(pow(a[i],2)+2.0*r[i]*delta[i]+pow(delta[i],2))/(2.0*a[i]*r[i]+2.0*a[i]*delta[i]));
			Dlg.spec_gasatt(f,Th, rou_h,  Ph, &gama);
			TotalAtt = a[i]*gama+TotalAtt;
		}
	}
	*attgas = TotalAtt;
}

 
////////////////////////////////中间程序////////////////////////////////////
//计算大气衰减率
void CMultDlg::spec_gasatt(double freq, double T, double Rou, double P, double *gama)
{
	//函数功能:计算得到给定温湿压条件下的特征衰减;
	//输入参数:频率（GHz），温度（K），水汽密度（读数据库g/m^3）,压强（Pa）
    //输出参数:特征衰减; 	
	double p;               //干空气压强(总压强为水蒸气压强和干空气压强之和即P=e+p);
	double e;               //水蒸气压强(随高度的变化会输入)
//	double T;              //温度(k);
    double cita;            //温度(转换);
//	double Rou;             //水蒸气密度;
    double SO;          //氧气吸收线强度;
	double FO;          //氧气吸收线线性因子;
    double SW;          //水蒸气吸收线强度;
    double FW;          //水蒸气吸收线线性因子;
	double TO_att=0.0;      //氧气衰减;
    cita=300/T;
//水蒸汽压强的计算；	
	e=Rou*T/216.7;
	p=P-e; //干空气压强；
//定义变量
	double dtfo;            //氧气线宽;
    double dtfw;            //水蒸气线宽;
	double dalto;           //氧气修正因子;
    double daltw;           //水气修正因子;
	double NDf;             //大气压力造成的氮气吸收连续带; 
	double NF;              //复折射率的虚部;
    
	//氧气衰减谱数据
	double f0[44]={50.474238,50.987749,51.503350,52.02141,52.542394,53.066907,53.595749,54.13,54.671159,55.221367,55.783802,56.264775,56.363389,56.968206,57.612484,58.323877,58.44659,59.164207,59.590983,60.306061,60.434776,61.15056,
		61.800154,62.411215,62.48626,62.997977,63.568518,64.127767,64.678903,65.224071,65.764772,66.302091,66.83683,67.369598,67.900867,68.431005,68.960311,118.750343,368.498350,424.763124,487.249370,715.393150,773.839675,834.14533};
	double a1[44]={0.94,2.46,6.08,14.14,31.02,64.10,124.7,228.00,391.80,631.6,953.5,548.9,1344.00,1763.00,2141.00,2386.0,1457.0,2404.0,2112.0,2124.0,2461.0,2504.0,2298.0,1933.0,1517.0,1503.0,1087.0,733.0,463.0,274.0,153.0,80.09,39.46,18.32,8.01,3.30,1.28,945.0,67.9,638.0,235.0,99.60,671.0,180.0}; 
	double a2[44]={9.694,8.694,7.744,6.844,6.004,5.224,4.484,3.814,3.194,2.624,2.119,0.015,1.66,1.26,0.915,0.626,0.084,0.391,0.212,0.212,0.391,0.626,0.915,1.26,0.083,1.665,2.115,2.62,3.195,3.815,4.485,5.225,6.005,6.845,7.745,8.695,9.695,0.009,0.049,0.044,0.049,0.145,0.13,0.147}; 
	double a3[44]={8.9,9.1,9.4,9.7,9.9,10.2,10.5,10.7,11,11.3,11.7,17.3,12.0,12.4,12.8,13.3,15.2,13.9,14.3,14.5,13.6,13.1,12.7,12.3,15.4,12.0,11.7,11.3,11.0,10.7,10.5,10.2,9.9,9.7,9.4,9.2,9.0,16.3,19.2,19.3,19.2,18.1,18.2,18.1};
	double a4[44]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.6,0.6,0.6,0.6,0.6,0.6};
	double a5[44]={2.4,2.2,1.97,1.66,1.36,1.31,2.3,3.35,3.74,2.58,-1.66,3.90,-2.97,-4.16,-6.13,-2.05,7.48,-7.22,7.65,-7.05,6.97,1.04,5.7,3.6,-4.98,2.39,1.08,-3.11,-4.21,-3.75,-2.67,-1.68,-1.69,-2.0,-2.28,-2.40,-2.5,-0.36,0,0,0,0,0,0};
	double a6[44]={7.9,7.8,7.74,7.64,7.51,7.14,5.84,4.31,3.05,3.39,7.05,-1.13,7.53,7.42,6.97,0.51,-1.46,2.66,-0.9,0.81,-3.24,-0.67,-7.61,-7.77,0.97,-7.68,-7.06,-3.32,-2.98,-4.23,-5.75,-7.0,-7.35,-7.44,-7.53,-7.6,-7.65,0.09,0,0,0,0,0,0};
    //水气衰减谱数据
	double fw[35]={22.23508,67.80396,119.99594,183.310091,321.225644,325.152919,336.222601,380.197372,390.134508,437.346667,439.150812,443.018295,448.001075,470.888947,474.689127,488.491133,503.569532,504.482692,547.67644,552.02096,556.936002,620.700807,
		645.866155,658.00528,752.033227,841.053973,859.962313,899.306675,902.616173,906.207325,916.171582,923.118427,970.315022,987.926764,1780.0};
	double b1[35]={0.113,0.0012,0.0008,2.42,0.0483,1.499,0.0011,11.52,0.0046,0.065,0.9218,0.1976,10.32,0.3297,1.262,0.252,0.039,0.013,9.701,14.77,487.4,5.012,0.0713,0.3022,239.6,0.014,0.1472,0.0605,0.0426,0.1876,8.34,0.0869,8.972,132.1,22300.0}; 
	double b2[35]={2.143,8.735,8.356,0.668,6.181,1.54,9.829,1.048,7.35,5.05,3.596,5.05,1.405,3.599,2.381,2.853,6.733,6.733,0.114,0.114,0.159,2.20,8.58,7.82,0.396,8.18,7.989,7.917,8.432,5.111,1.442,10.22,1.92,0.258,0.952}; 
	double b3[35]={28.11,28.58,29.48,30.50,23.03,27.83,26.93,28.73,21.52,18.45,21.00,18.60,26.32,21.52,23.55,26.02,16.12,16.12,26.00,26.00,32.10,24.38,18.00,32.10,30.6,15.9,30.6,29.85,28.65,24.08,26.7,29.0,25.5,29.85,176.2};
	double b4[35]={0.69,0.69,0.7,0.64,0.67,0.68,0.69,0.54,0.63,0.6,0.63,0.6,0.66,0.66,0.65,0.69,0.61,0.61,0.7,0.7,0.69,0.71,0.6,0.69,0.68,0.33,0.68,0.68,0.7,0.7,0.7,0.7,0.64,0.68,0.5};
	double b5[35]={4.8,4.93,4.78,5.3,4.69,4.85,4.74,5.38,4.81,4.23,4.29,4.23,4.84,4.57,4.65,5.04,3.89,4.01,4.5,4.5,4.11,4.68,4.0,4.14,4.09,5.76,4.09,4.53,5.1,4.7,4.78,5.0,4.94,4.55,30.5};
	double b6[35]={1.0,0.82,0.79,0.85,0.54,0.74,0.61,0.89,0.55,0.48,0.52,0.5,0.67,0.65,0.64,0.72,0.43,0.45,1.0,1.0,1.0,0.68,0.5,1.0,0.84,0.45,0.84,0.9,0.95,0.53,0.78,0.8,0.67,0.9,1.0};
	double A1,A2,B1,B2;
	double NFo=0;             //氧气SF乘积的和;
	double NFO;

	//计算氧气衰减效应	
	//频率小于指定值时，考虑所有氧气谱线
	if (freq<=118.750343)
	{
		for (int i=0;i<44;i++)
		{
			dtfo=a3[i]*1E-4*(p*pow(cita,0.8-a4[i])+1.1*e*cita);//线宽
			dtfo=sqrt(pow(dtfo,2)+2.25E-6); //修正
			dalto=(a5[i]+a6[i]*cita)*1E-4*P*pow(cita,0.8);//修正因子
			A1=(dtfo-dalto*(f0[i]-freq))/(pow(f0[i]-freq,2)+pow(dtfo,2));
			A2=(dtfo-dalto*(f0[i]+freq))/(pow(f0[i]+freq,2)+pow(dtfo,2));
			FO=freq/f0[i]*(A1+A2);
			SO=a1[i]*pow(10,-7)*p*pow(cita,3)*exp(a2[i]*(1-cita)); //线强度
			NFO=SO*FO;
			NFo=NFo+NFO;	
		}
	}
	else if (freq>118.750343)//高于指定频率，只考虑60GHz以上谱线
	{
		for (int i=19;i<44;i++)
		{
			dtfo=a3[i]*1E-4*(p*pow(cita,0.8-a4[i])+1.1*e*cita);//线宽
			dtfo=sqrt(pow(dtfo,2)+2.25E-6); //修正
			dalto=(a5[i]+a6[i]*cita)*1E-4*P*pow(cita,0.8);//修正因子
			A1=(dtfo-dalto*(f0[i]-freq))/(pow(f0[i]-freq,2)+pow(dtfo,2));
			A2=(dtfo-dalto*(f0[i]+freq))/(pow(f0[i]+freq,2)+pow(dtfo,2));
			FO=freq/f0[i]*(A1+A2);
			SO=a1[i]*pow(10,-7)*p*pow(cita,3)*exp(a2[i]*(1-cita)); //线强度
			NFO=SO*FO;
			NFo=NFo+NFO;			
		}
	}	
//计算水蒸气衰减效应;
	double NFw=0;             //水蒸气SF乘积的和;
	double NFW;
	double dtfw1;       //中间变量;
	for (int j=0;j<35;j++)
	{		
		dtfw1=b3[j]*1E-4*(p*pow(cita,b4[j])+b5[j]*e*pow(cita,b6[j])); //线宽
		dtfw=0.535*dtfw1+sqrt(0.217*pow(dtfw1,2)+2.1316E-12*pow(fw[j],2)/cita);	//修正后的线宽		
		daltw=0; //修正因子;		
		B1=(dtfw-daltw*(fw[j]-freq))/(pow(fw[j]-freq,2)+pow(dtfw,2));
		B2=(dtfw-daltw*(fw[j]+freq))/(pow(fw[j]+freq,2)+pow(dtfw,2));

		FW=freq/fw[j]*(B1+B2);
        SW=b1[j]*1E-1*e*pow(cita,3.5)*exp(b2[j]*(1-cita));
		NFW=SW*FW;
        NFw=NFw+NFW;		
	}
    
	//大气压力造成的大气吸收和deby频谱的干燥吸收连续带;
	double C1,C2,d;
	d=5.6E-4*P*pow(cita,0.8);
    C1=6.14E-5/(d*(1+pow(freq/d,2)));
    C2=1.4E-12*p*pow(cita,1.5)/(1+1.9E-5*pow(freq,1.5));
    NDf=freq*p*pow(cita,2)*(C1+C2);
	
	//求解频率依赖的复折射率的虚部;
    NF=NFo+NFw+NDf;
	
	//大气特征衰减;
    *gama=0.1820*freq*NF;
}


///////////////////////////云雾衰减计算函数CalAttCloud////////////////////////////////////////
void CMultDlg::CalAttCloud(double Ts, double Ps, double rou0, double rain_h,double plat_h1,double plat_h2,double f,double theta0,double LW,double *attcloud)
{
		
		//输入参数：雨顶高度rain_h(m),发射高度 plat_h1（m），接收高度 plat_h2（m），温度（K,计算云衰减，此处为273.15K），频率（GHz）,仰角，液态水含量(kg/m^2)基于概率读数据库
	//输出，云衰减
	//模型常系数
	CMultDlg Dlg;

	double pi=3.1415926;
	
	double T = 273.15;//设定为273.15，为计算云衰减，当计算雾衰减时应采用实际输入值
	rain_h=rain_h/1000;
	plat_h1=plat_h1/1000;
	plat_h2=plat_h2/1000;

	
	double Thh=Ts-6.5*rain_h;
	double Ph=Ps*pow(Ts/(Ts-6.5*rain_h),-34.163/6.5);
	double rou_h=rou0*exp(-rain_h*0.5);
	double e = rou_h*Thh/216.7;
	double Nh = 77.6/Thh*(Ph+4810.0*e/Thh);
	double nh = 1+Nh*1E-6;
	
	double Ths=Ts-6.5*plat_h1;
	double Phs=Ps*pow(Ts/(Ts-6.5*plat_h1),-34.163/6.5);
	double rou_hs=rou0*exp(-plat_h1*0.5);
	double es = rou_hs*Ths/216.7;
	double Nhs = 77.6/Ths*(Phs+4810.0*es/Ths);
	double nhs = 1+Nhs*1E-6;

	
	double costh2=nhs*(6371000+plat_h1*1000)*cos(theta0*pi/180)/(nh*(6371000+plat_h2*1000));
	double thetah=acos(costh2)*180/pi;
	

	double cloud_h = 6;//定义云顶高度为6km
	double Th=300/T;
	double e0=77.6+103.3*(Th-1);
	double e1=5.48;
	double e2=3.51;
	double Re = 8500;
////////////////////
	double fp;
	double fs;
	double fpp;
	double fss;
	double er;//复介电常数的实部
	double ei;//复介电常数的虚部
	double yeta;
	double Kl;
	double eff_h;
	fp=20.09-142.4*(Th-1)+294*(Th-1)*(Th-1);
	fs=590-1500*(Th-1);
	fpp = pow((f/fp),2);
	fss = pow((f/fs),2);
	er = (e0-e1)/(1+fpp)+(e1-e2)/(1+fss)+e2;
	ei = f*(e0-e1)/(fp*(1+fpp))+f*(e1-e2)/(fs*(1+fss));
	yeta = (2+er)/ei;
	Kl = 0.819*f/(ei*(1+pow(yeta,2)));	
    eff_h = cloud_h-rain_h;
	if(plat_h2>plat_h1)
	{
		if (plat_h1>cloud_h && plat_h2>cloud_h)//收发皆在云顶高度以上
			LW = LW*0;
		else if(plat_h1<rain_h && plat_h2<rain_h)//收发皆在雨顶高度以下
			LW = LW*0;
		else if(plat_h1<rain_h && rain_h<plat_h2 && plat_h2<cloud_h)//低高度平台在雨顶高度以下，高平台在雨顶高度与云顶高度之间
		{
			LW = (plat_h2-rain_h)/(cloud_h-rain_h)*LW;
		    eff_h = plat_h2-rain_h;
		}
		else if(rain_h<=plat_h1 && plat_h1<=cloud_h && rain_h<=plat_h2 && plat_h2<=cloud_h)//收发皆在雨顶高度与云顶高度之间
		{
			LW = (plat_h2-plat_h1)/(cloud_h-rain_h)*LW;
            eff_h = plat_h2-plat_h1;
		}
		else if(rain_h<plat_h1 && plat_h1<cloud_h && plat_h2>cloud_h)//低平台在雨顶高度与云顶高度之间，高平台在云顶高度以上
		{
			LW = (cloud_h-plat_h1)/(cloud_h-rain_h)*LW;
            eff_h = (cloud_h-plat_h1);
		}
	}
	else if((plat_h2<plat_h1))
	{
		if (plat_h1>cloud_h && plat_h2>cloud_h)//收发皆在云顶高度以上
			LW = LW*0;
		else if(plat_h1<rain_h && plat_h2<rain_h)//收发皆在雨顶高度以下
			LW = LW*0;
		else if(plat_h2<rain_h && rain_h<plat_h1 && plat_h1<cloud_h)//低高度平台在雨顶高度以下，高平台在雨顶高度与云顶高度之间
		{
			LW = (plat_h1-rain_h)/(cloud_h-rain_h)*LW;
            eff_h = (plat_h1-rain_h);
		}
		else if(rain_h<=plat_h2 && plat_h2<=cloud_h && rain_h<=plat_h1 && plat_h1<=cloud_h)//收发皆在雨顶高度与云顶高度之间
		{
			LW = (plat_h1-plat_h2)/(cloud_h-rain_h)*LW;
		    eff_h = (plat_h1-plat_h2);
		}
		else if(rain_h<plat_h2 && plat_h2<cloud_h && plat_h1>cloud_h)//低平台在雨顶高度与云顶高度之间，高平台在云顶高度以上
		{
			LW = (cloud_h-plat_h2)/(cloud_h-rain_h)*LW;
		    eff_h = (cloud_h-plat_h2);
		}
	}

	if (theta0>=5)
	{
		if(plat_h1<rain_h)
			*attcloud=Kl*LW/sin(thetah*pi/180);
		else
			*attcloud=Kl*LW/sin(theta0*pi/180);
	}
	else
	{
		if(plat_h1<rain_h)
		{
			LW = 2*LW/(pow(sin(thetah*pi/180)*sin(thetah*pi/180)+2*eff_h/Re,0.5)+sin(thetah*pi/180));	
			*attcloud=Kl*LW;	
		}
		else
		{
			LW = 2*LW/(pow(sin(theta0*pi/180)*sin(theta0*pi/180)+2*eff_h/Re,0.5)+sin(theta0*pi/180));	
			*attcloud=Kl*LW;	
		}
	}


//	if (theta0>=5)
//	{
//		*attcloud=Kl*LW/sin(theta0*pi/180);
//	}
//	else
//	{
//		LW = 2*LW/(pow(sin(theta0*pi/180)*sin(theta0*pi/180)+2*eff_h/Re,0.5)+sin(theta0*pi/180));	
//		*attcloud=Kl*LW;
//	}
}

///////////////////////////降雨衰减计算函数CalAttRain////////////////////////////////////////

void CMultDlg::CalAttRain(double lats, double theta0, double ipol, double f, double hs, double hr,double ha, double R001, double p, double *attrain)
{
/*	CMultDlg Dlg;
	double pi=3.1415926;
	
	p=100-p;
	hs=hs/1000;
	ha=ha/1000;
	hr=hr/1000;
	if(p>5)
		*attrain=0;
	else
	{
		//由极化方式确定极化角
	double tal;
	if(ipol==0)
		tal = 0;
	else if(ipol==1)
		tal = 90;
	else if (ipol==2)
		tal = 45;
//角度转换
	//tal = tal*pi/180;
//	theta0 = theta0*pi/180;
	

//特征衰减
	double k,arfa;
	double a_kh[4] = {-5.3398,-0.35351,-0.23789,-0.94158};
	double a_kv[4] = {-3.80595,-3.44965,-0.39902,0.50167};
	double b_kh[4] = {-0.10008,1.2697,0.86036,0.64552};
	double b_kv[4] = {0.56934,-0.22911,0.73042,1.07319};
	double c_kh[4] = {1.13098,0.454,0.15354,0.16817};
	double c_kv[4] = {0.81061,0.51059,0.11899,0.27195};
	double a_ah[5] = {-0.14318,0.29591,0.32177,-5.3761,16.1721};
	double a_av[5] = {-0.07771, 0.56727, -0.20238, -48.2991, 48.5833};
	double b_ah[5] = {1.82442, 0.77564, 0.63773, -0.9623, -3.2998};
	double b_av[5] = {2.3384, 0.95545, 1.1452, 0.791669, 0.791459};
	double c_ah[5] = {-0.55187, 0.19822, 0.13164, 1.47828, 3.4399};
	double c_av[5] = {-0.76284, 0.54039, 0.26809, 0.116226, 0.116479};
	double ma_ah = 0.67849;
	double ca_ah = -1.95537;
	double ma_av = -0.053739;
	double ca_av = 0.83433;
	double mk_kh = -0.18961;
	double ck_kh = 0.71147;
	double mk_kv = -0.16398;
	double ck_kv = 0.63297;
	double sum_kh=0,sum_kv=0,sum_ah=0,sum_av=0;
	double kH,kV,aH,aV;
	for(int jh=0;jh<=3;jh++)
	{
		sum_kh = sum_kh+a_kh[jh]*exp(-pow((log10(f)-b_kh[jh])/c_kh[jh],2.0));
        sum_kv = sum_kv+a_kv[jh]*exp(-pow((log10(f)-b_kv[jh])/c_kv[jh],2.0));
	}
	sum_kh = sum_kh+mk_kh*log10(f)+ck_kh;
	sum_kv = sum_kv+mk_kv*log10(f)+ck_kv;
	kH = pow(10,sum_kh);
	kV = pow(10,sum_kv);
	for(int j=0;j<=4;j++)
	{ 
		sum_ah = sum_ah+a_ah[j]*exp(-pow((log10(f)-b_ah[j])/c_ah[j],2.0));
		sum_av = sum_av+a_av[j]*exp(-pow((log10(f)-b_av[j])/c_av[j],2.0));
	}
	aH = sum_ah+ma_ah*log10(f)+ca_ah;
	aV = sum_av+ma_av*log10(f)+ca_av;
	k = (kH+kV+(kH-kV)*pow(cos(theta0*pi/180.0),2.0)*cos(2*tal*pi/180.0))/2.0;
	arfa = (kH*aH+kV*aV+(kH*aH-kV*aV)*pow(cos(theta0*pi/180.0),2.0)*cos(2.0*tal*pi/180.0))/(2.0*k);
	double rR = k*pow(R001,arfa);
//雨衰减
// 比较雨顶高度与飞行器高度
	if(hr>=ha)
	{
		hr = ha;
	}

//斜路径长度
	double LS,LG;
	if(theta0>=5)
		LS = (hr-hs)/sin(theta0*pi/180.0);
	else
		LS = 2*(hr-hs)/(sqrt(pow(sin(theta0*pi/180.0),2)+2*(hr-hs)/8500)+sin(theta0*pi/180.0));
//斜路径水平投影	
	LG = LS*cos(theta0*pi/180.0);
	double rp = 3.78*pow(R001,(-0.56+1.51/LS))*(1-0.85*(pow(p,0.065))/(1+0.12*LS));
if(rp>2.5)
    rp = 2.5;


double A001 = k*pow((rp*R001),arfa)*LS;


	double bita;
	if((p>=1)||(lats>=36))
		bita=0;
	else
	{
		if((p<1)&&(lats<36)&&(theta0>=25))
			bita = -0.005*(lats-36.0);
		else
			bita = -0.005*(lats-36.0)+1.8-4.25*sin(theta0*pi/180.0);
	}
	*attrain = A001*pow(p/0.01,-(0.655+0.033*log(p)-0.045*log(A001)-bita*(1.0-p)*sin(theta0*pi/180.0)));
	}*/
	CMultDlg Dlg;
	double pi=3.1415926;
	
	p=100-p;
	hs=hs/1000;
	ha=ha/1000;
	hr=hr/1000;
	if(p>5)
		*attrain=0;
	else
	{
		//由极化方式确定极化角
	double tal;
	if(ipol==0)
		tal = 0;
	else if(ipol==1)
		tal = 90;
	else if (ipol==2)
		tal = 45;
//角度转换
	//tal = tal*pi/180;
//	theta0 = theta0*pi/180;
	

//特征衰减
	double k,arfa;
	double a_kh[4] = {-5.3398,-0.35351,-0.23789,-0.94158};
	double a_kv[4] = {-3.80595,-3.44965,-0.39902,0.50167};
	double b_kh[4] = {-0.10008,1.2697,0.86036,0.64552};
	double b_kv[4] = {0.56934,-0.22911,0.73042,1.07319};
	double c_kh[4] = {1.13098,0.454,0.15354,0.16817};
	double c_kv[4] = {0.81061,0.51059,0.11899,0.27195};
	double a_ah[5] = {-0.14318,0.29591,0.32177,-5.3761,16.1721};
	double a_av[5] = {-0.07771, 0.56727, -0.20238, -48.2991, 48.5833};
	double b_ah[5] = {1.82442, 0.77564, 0.63773, -0.9623, -3.2998};
	double b_av[5] = {2.3384, 0.95545, 1.1452, 0.791669, 0.791459};
	double c_ah[5] = {-0.55187, 0.19822, 0.13164, 1.47828, 3.4399};
	double c_av[5] = {-0.76284, 0.54039, 0.26809, 0.116226, 0.116479};
	double ma_ah = 0.67849;
	double ca_ah = -1.95537;
	double ma_av = -0.053739;
	double ca_av = 0.83433;
	double mk_kh = -0.18961;
	double ck_kh = 0.71147;
	double mk_kv = -0.16398;
	double ck_kv = 0.63297;
	double sum_kh=0,sum_kv=0,sum_ah=0,sum_av=0;
	double kH,kV,aH,aV;
	for(int jh=0;jh<=3;jh++)
	{
		sum_kh = sum_kh+a_kh[jh]*exp(-pow((log10(f)-b_kh[jh])/c_kh[jh],2.0));
        sum_kv = sum_kv+a_kv[jh]*exp(-pow((log10(f)-b_kv[jh])/c_kv[jh],2.0));
	}
	sum_kh = sum_kh+mk_kh*log10(f)+ck_kh;
	sum_kv = sum_kv+mk_kv*log10(f)+ck_kv;
	kH = pow(10,sum_kh);
	kV = pow(10,sum_kv);
	for(int j=0;j<=4;j++)
	{ 
		sum_ah = sum_ah+a_ah[j]*exp(-pow((log10(f)-b_ah[j])/c_ah[j],2.0));
		sum_av = sum_av+a_av[j]*exp(-pow((log10(f)-b_av[j])/c_av[j],2.0));
	}
	aH = sum_ah+ma_ah*log10(f)+ca_ah;
	aV = sum_av+ma_av*log10(f)+ca_av;
	k = (kH+kV+(kH-kV)*pow(cos(theta0*pi/180.0),2.0)*cos(2*tal*pi/180.0))/2.0;
	arfa = (kH*aH+kV*aV+(kH*aH-kV*aV)*pow(cos(theta0*pi/180.0),2.0)*cos(2.0*tal*pi/180.0))/(2.0*k);
	double rR = k*pow(R001,arfa);
//雨衰减
// 比较雨顶高度与飞行器高度
	if(hr>=ha)
	{
		hr = ha;
	}

//斜路径长度
	double LS,LG;
	if(theta0>=5)
		LS = (hr-hs)/sin(theta0*pi/180.0);
	else
		LS = 2*(hr-hs)/(sqrt(pow(sin(theta0*pi/180.0),2)+2*(hr-hs)/8500)+sin(theta0*pi/180.0));
//斜路径水平投影	
	LG = LS*cos(theta0*pi/180.0);
//水平缩减因子
	double r001 = 1/(1+0.78*sqrt(LG*rR/f)-0.38*(1-exp(-2*LG)));
//垂直调整因子
	double x;
	if(lats<36)		
		x = 36-lats;
	else		
		x = 0;
	double kesi = atan((hr-hs)/(LG*r001))*180.0/pi;//计算并转换为角度
	double LR;
	if(kesi>theta0)
		LR = LG*r001/cos(theta0*pi/180.0);
	else   
		LR = (hr-hs)/sin(theta0*pi/180.0);
	double v001 = 1.0/(1.0+sqrt(sin(theta0*pi/180.0))*(31.0*(1.0-exp(-theta0/(1.0+x)))*sqrt(LR*rR)/pow(f,2.0)-0.45));
//有效路径长度
	double LE = LR*v001;
//0.01%的雨衰减
	double A001 = rR*LE;
//p%的雨衰减
	double bita;
	if((p>=1)||(lats>=36))
		bita=0;
	else
	{
		if((p<1)&&(lats<36)&&(theta0>=25))
			bita = -0.005*(lats-36.0);
		else
			bita = -0.005*(lats-36.0)+1.8-4.25*sin(theta0*pi/180.0);
	}
	*attrain = A001*pow(p/0.01,-(0.655+0.033*log(p)-0.045*log(A001)-bita*(1.0-p)*sin(theta0*pi/180.0)));
	}
}
                       

//////////////多径与闪烁衰落函数CalAttMS的中间函数SintAtt////////////////////////////////////////
void CMultDlg::SintAtt(int season, double Ts, double rou0, double lat, double theta, double f, double p0, double d, double yita, double *attsint, double *newp)
{
//输入参数：季节参数，地面节点纬度，路径仰角、频率、时间概率、天线直径、天线效率
//输出参数：闪烁衰落和多径衰落
	//常数
	
	double e=rou0*Ts/216.7;//折射率湿项，N 
	double nwet=3.732*pow(10,5)*e/(Ts*Ts);
    double pi=3.1415926;

/////////////////////////////////////////////
	double Gw=0;
    double Kw=0;
	double dG=0;
	double pL=0;
	double Clat=0;
    double C0=0;
	double dtref=0;
	double L=0;
	double Deff=0;
	double x=0;
	double gx=0;
	double dt=0;
	double ap=0;
	double Abs=0;
	double At=25;
	double A63=0;
	double pt=0;
	double qq=0;
	double s0=0;
	double qt=0;
	double q=0;
	double p=0;


	if(season==1)
		pL=20;
	if(season==2)
		pL=5;
	if(season==3)
		pL=5;
	if(season==4)
		pL=1;


	if(theta<5)
	{
		A63=2.27-1.16*log10(1+theta*1000*pi/180);

		if(lat<=53)
			Clat=0;
		else if(lat>53 && lat<60)
			Clat=-53+lat;
		else
			Clat=7;

		C0=76;
		Kw=pow(pL,1.5)*pow(10,0.1*(C0+Clat));
		if (lat<=45)
			dG=-1.8-5.6*log10(1.1+pow(fabs(cos(2*lat*pi/180)),0.7))+4.5*log10(1+theta*pi*1000/180);
		else
            dG=-1.8-5.6*log10(1.1-pow(fabs(cos(2*lat*pi/180)),0.7))+4.5*log10(1+theta*pi*1000/180);
		
		Gw=10*log10(Kw)-92;
		Gw=Gw-dG;
		*attsint=Gw+92+9*log10(f)-55*log10(1+theta*pi*1000/180)-10*log10(p0);
		*newp=p0;

		if(*attsint<A63+At)
		{
			pt=Kw*pow(10,-0.1*dG)*pow(f,0.9)*pow(1+theta*1000*pi/180,-5.5)*pow(10,-At/10);
			p=pow(10,-0.1*A63+log10(pt));
			qq=-20*log10(-log((100-p)/100))/At;
			s0=-1.6-3.2*log10(f)+4.2*log10(1+theta*1000*pi/180);
			qt=(qq-2)/((1+0.3*pow(10,-At/20))*pow(10,-0.016*At))-s0*(pow(10,-At/20)+At/800);
			if(qt<0)
			{
				At=35;
				pt=Kw*pow(f,0.9)*pow(1+theta*1000*pi/180,-5.5)*pow(10,-At/10);
				p=pow(10,-0.1*A63+log10(pt));
				qq=-20*log10(-log((100-p)/100))/At;
				s0=-1.6-3.2*log10(f)+4.2*log10(1+theta*1000*pi/180);
				qt=(qq-2)/((1+0.3*pow(10,-At/20))*pow(10,-0.016*At))-s0*(pow(10,-At/20)+At/800);
			}
			if(*attsint>A63)
			{
				q=2+pow(10,-0.016*(*attsint-A63))*(1+0.3*pow(10,-(*attsint-A63)/20))*(qt+s0*(pow(10,-(*attsint-A63)/20)+(*attsint-A63)/800));
				*newp=100*(1-exp(-pow(10,-q*(*attsint-A63)/20)));
			}
			else if(*attsint<=A63)
			{
				double Eref=0;
				Eref=A63-*attsint;
				double A001=0;
				A001=Gw+92+9*log10(f)-55*log10(1+theta*pi*1000/180)-10*log10(0.01);
				
				if(Eref>10)
				{
					*newp=pow(10,((-1.7+0.2*A001-Eref)/3.5));
				}
				else if(Eref>0 && Eref<=10)
				{
					double pw1=0;
					double qe1=0;
					double qs=0;
					double qe=0;
					pw1=100-pow(10,(-1.7+0.2*A001-10)/3.5);
					qe1=-20*(log10(-log(1-(100-pw1)/58.21)))/10;
					qs=2.05*qe1-20.3;
					qe=8+(1+0.3*pow(10,-Eref/20))*(pow(10,-0.7*Eref/20))*(qs+12*(pow(10,-Eref/20)+Eref/800));
					*newp=58.21*(1-exp(-pow(10,-qe*Eref/20)));
				}
				*attsint=-Eref;
			}
		}
	}
	else
	{
		dtref=3.6*pow(10,-3)+pow(10,-4)*nwet;
        //有效路径长度
		L=2*1000/(sqrt(pow(sin(theta*pi/180.0),2.0)+2.35*pow(10,-4))+sin(theta*pi/180.0));
        //有效天线直径
		Deff=sqrt(yita)*d;
        //天线平均因子
	    x=1.22*pow(Deff,2.0)*f/L;
	    gx=sqrt(3.86*pow(x*x+1.0,11.0/12.0)*sin(atan(1.0/x)*11.0/6.0)-7.08*pow(x,5.0/6.0));
        //信号标准偏差
	    dt=dtref*pow(f,7.0/12.0)*gx/pow(sin(theta*pi/180),1.2);
        //时间概率因子
		ap=-0.061*pow(log10(p0),3.0)+0.072*pow(log10(p0),2.0)-1.71*log10(p0)+3.0;
        //闪烁衰落
		*attsint=ap*dt;
		*newp=p0;
	}
}

/////////多径小尺度衰落与多径信道模型函数CalAttMult的中间函数ontrt，ocdiv////////////////////////////////////////
void CMultDlg::ontrt(dxcomplex *a, int n)
{//复数的N次方根
    double r,q,t;
    if (n<1) return;
    q=atan2(a->imag,a->real);
    r=sqrt(a->real*a->real+a->imag*a->imag);
    if (r+1.0!=1.0)
	{ r=(1.0/n)*log(r); r=exp(r);}
	
    t=q/n;
    a->real=r*cos(t); 
	a->imag=r*sin(t);
    return;
}

void CMultDlg::ocdiv(dxcomplex *a, dxcomplex *b)
{//复数除法
	double p,q,s,w;
    p=a->real*b->real; q=-a->imag*b->imag; s=(a->real+a->imag)*(b->real-b->imag);
    w=b->real*b->real+b->imag*b->imag;
    if (w+1.0==1.0) 
	{ a->real=1.0e+35*a->real/fabs(a->real);
	a->imag=1.0e+35*a->imag/fabs(a->imag);
	}
    else
	{ a->real=(p-q)/w; a->imag=(s-p-q)/w; }
    return;
}

////////////////////多径与闪烁衰落函数CalAttMS/////////////////////////////////////////////////////////////////////////////////
void CMultDlg::CalAttMS(int season, double Ts, double rou0, double lats, double theta, double f, double pout, double d, double yita, double *AMSout)
{
	CMultDlg Dlg;

	double attsint=0;
	double newp=0;
	double p0=0;

	int npp=40;
	double* pinp = NULL;
	pinp = new double[npp];
	memset(pinp, 0, sizeof(double) * npp); 

	double* attMS = NULL;
	attMS = new double[npp];
	memset(attMS, 0, sizeof(double) * npp); 


	pinp[0]=99.999;
	pinp[1]=99.998;
	pinp[2]=99.997;
	pinp[3]=99.995;
	pinp[4]=99.99;
	pinp[5]=99.98;
	pinp[6]=99.97;
	pinp[7]=99.95;
	pinp[8]=99.9;
	pinp[9]=99.8;
	pinp[10]=99.7;
	pinp[11]=99.5;
	pinp[12]=99;
	pinp[13]=98;
	pinp[14]=97;
	pinp[15]=95;
	pinp[16]=90;
	pinp[17]=80;
	pinp[18]=70;
	pinp[19]=60;
	pinp[20]=50;
	pinp[21]=30;
	pinp[22]=20;
	pinp[23]=10;
	pinp[24]=5;
	pinp[25]=3;
	pinp[26]=2;
	pinp[27]=1;
	pinp[28]=0.5;
	pinp[29]=0.3;
	pinp[30]=0.2;
	pinp[31]=0.1;
	pinp[32]=0.05;
	pinp[33]=0.03;
	pinp[34]=0.02;
	pinp[35]=0.01;
	pinp[36]=0.005;
	pinp[37]=0.003;
	pinp[38]=0.002;
	pinp[39]=0.001;

    bool ifind=false;

	for(int i=0;i<npp;i++)
	{
		p0=pinp[i];
		Dlg.SintAtt(season, Ts, rou0, lats, theta, f, p0, d, yita, &attsint, &newp);
		attMS[i]=attsint;
		pinp[i]=newp;
		if(100-pout==pinp[i])
		{
			*AMSout=attMS[i];
			ifind =true;
			break;
		}
	}
	for(i=1;i<npp;i++)
	{
		
		if(100-pinp[i-1]<pout && 100-pinp[i]>pout)
		{
			*AMSout=(attMS[i]-attMS[i-1])*(pout-100+pinp[i-1])/(pinp[i-1]-pinp[i])+attMS[i-1];
			ifind=true;
			break;
		}

	}
	if(ifind==false)
	{
		double dp=fabs(100-pinp[0]-pout);
		*AMSout=attMS[0];
		for(i=1;i<npp;i++)
		{
			if(fabs(100-pinp[i]-pout)<dp)
			{
				dp=fabs(100-pinp[i]-pout);
				*AMSout=attMS[i];
			}
		}
		ifind=true;
	}
}
////////////////////新的多径与闪烁衰落函数newCalAttMS/////////////////////////////////////////////////////////////////////////////////
void CMultDlg::newCalAttMS(double rou0, double Ts, int season, double hs, double lat, double theta, double f, double p0, double d, double yita, double *attsint)
{
//输入参数：季节参数，地面节点纬度，路径仰角、频率、时间概率、天线直径、天线效率
//输出参数：闪烁衰落和多径衰落
	//常数
	p0=100-p0;
	if(p0>50)
		*attsint=0;
	else
	{
		double e=rou0*Ts/216.7;
		double nwet=3.732*pow(10,5)*e/(Ts*Ts);//折射率湿项，N 
	////	nwet=57.4;/////////////////////////////////////////////////////////////////
		double pi=3.1415926;

/////////////////////////////////////////////
		double Gw=0;
		double Kw=0;
		double dG=0;
		double pL=0;
		double Clat=0;
		double C0=0;
		double dtref=0;
		double L=0;
		double Deff=0;
		double x=0;
		double gx=0;
		double dt=0;
		double ap=0;
		double Abs=0;
		double A1=25;
		double A63=0;
		double pt=0;
		double qq=0;
		double s0=0;
		double qt=0;
		double q=0;
		double p=0;


		dtref=3.6*pow(10,-3)+pow(10,-4)*nwet;
		Deff=sqrt(yita)*d;




		if(season==1)
			pL=20;
		if(season==2)
			pL=5;
		if(season==3)
			pL=5;
		if(season==4)
			pL=1;

////////pL=9;///////////////////////////////////////
		if(theta<5)
		{

			if(lat<=53)
				Clat=0;
			else if(lat>53 && lat<60)
				Clat=-53+lat;
			else
				Clat=7;

			C0=76;/////////////+6*0.6///////////////////////////////////////////
			Kw=pow(pL,1.5)*pow(10,0.1*(C0+Clat));
			double v=0;
			if (lat<=45)
				v=-1.8-5.6*log10(1.1+pow(fabs(cos(2*lat*pi/180)),0.7));
			else
				v=-1.8-5.6*log10(1.1-pow(fabs(cos(2*lat*pi/180)),0.7));
		
			Gw=10*log10(Kw);
		
			*attsint=Gw-v+9*log10(f)-59*log10(1+theta*pi*1000/180)-10*log10(p0);

			if(*attsint<A1)
			{
				double theta1=0;
				double e=2.7183;
				theta1=pow(Kw*pow(10,-0.1*v)*pow(f,0.9)/(p0*pow(10,A1/10)),0.1681)-1;
				double A11=0;
				A11=-59.5*log10(e)/(1+theta1);
				double A2=0;
				double L5=2*1000/(sqrt(pow(sin(5*pi/180.0),2.0)+2.35*pow(10,-4))+sin(5*pi/180.0));
       			double x5=1.22*pow(Deff,2.0)*f/L5;
				double gx5=sqrt(3.86*pow(x5*x5+1.0,11.0/12.0)*sin(atan(1.0/x5)*11.0/6.0)-7.08*pow(x5,5.0/6.0));
				ap=-0.061*pow(log10(p0),3.0)+0.072*pow(log10(p0),2.0)-1.71*log10(p0)+3.0;
			
			
				double dt5=dtref*pow(f,7.0/12.0)*gx5/pow(sin(5*pi/180),1.2);
				A2=ap*dt5;
				double A21=0;
				double dgg=0;
				double dxth=0;
				double cc=0;
				cc=atan(1/x5)*11/6;
				dgg=(1770*(x5*x5+1)+2123*pow(x5,0.1667)*pow(x5*x5+1,0.9167)*(cos(cc)-x5*sin(cc)))/(12*pow(x5,0.1667)*(x5*x5+1)*(354*pow(x5,0.8333)-193*pow(x5*x5+1,0.9167)*sin(cc)));
				dxth=1.22*Deff*Deff*f*cos(5*pi/180)*(sin(5*pi/180)/sqrt(pow(sin(5*pi/180),2)+2.35*pow(10,-4))+1)/(2*1000);
				A21=A2*(dgg*dxth-1.2/tan(5*pi/180))/1000;
				double taos5=1/(1.728+0.5411*5+0.03723*5*5+hs/1000*(0.1815+0.06272*5+0.0138*5*5)+hs*hs*(0.01727+0.008288*5)/1000000);
				double theta2=(5+taos5)*1000*pi/180;
			
				double deta=theta2-theta1;
				double alpha=A11/A1;
				double beta=(log(A2/A1)-alpha*deta)/pow(deta,2);
				double gama=(A21-A2*(alpha+2*beta*deta))/(A2*deta*deta);
				*attsint=A1*exp(alpha*(theta*1000*pi/180-theta1)+beta*pow(theta*1000*pi/180-theta1,2)+gama*pow(theta*1000*pi/180-theta1,2)*(theta*1000*pi/180-theta2));
			}
		}
		else
		{
			//有效路径长度
			L=2*1000/(sqrt(pow(sin(theta*pi/180.0),2.0)+2.35*pow(10,-4))+sin(theta*pi/180.0));
			//天线平均因子
			x=1.22*pow(Deff,2.0)*f/L;
			gx=sqrt(3.86*pow(x*x+1.0,11.0/12.0)*sin(atan(1.0/x)*11.0/6.0)-7.08*pow(x,5.0/6.0));
			//信号标准偏差
			dt=dtref*pow(f,7.0/12.0)*gx/pow(sin(theta*pi/180),1.2);
			//时间概率因子
			ap=-0.061*pow(log10(p0),3.0)+0.072*pow(log10(p0),2.0)-1.71*log10(p0)+3.0;
			//闪烁衰落
			*attsint=ap*dt;
		}
	}
}



//////多径小尺度衰落与多径信道模型函数CalAttMult////////////////////////////////////////////////////////////////
void CMultDlg::CalAttMult(double theta, int season, double dsa, double hs, double ha, int ipol, double pout, double f, double* dt, double* damp, double* AMout,double* Lr)
{
	double pi = 3.1415926;
	double re=6371;//地球半径，km
	
	double wl = 0.3 / f;//波长
	double fko = 2 * pi / wl;//波数

	double re1=0;
	double re2=0;
	double grz=0;
	double rsa=0;
	double rsa1=0;
	double rsa2=0;
	double kre=1.3333*re;

	double epsilon=70;    //相对介电常数;
	double sigma=5;       //电导率


	double Y=0.25*kre*dsa*(hs-ha)/1000;
	double X=(kre*(hs+ha)/1000+dsa*dsa*0.25)/3;
	double fai=acos(Y/(X*sqrt(X)));
	double d31=2*sqrt(X)*cos((fai+4*pi)/3);
	double dt1=dsa/2-d31;
	double dr1=dsa/2+d31;
	double ht1=hs-dt1*dt1/(2000*kre);
	double hr1=ha-dr1*dr1/(2000*kre);
	double dr=2*ht1*hr1/(1000*dsa);
	double psi=(ht1+hr1)/dsa;

	double gamasp=0.0072*ha/1000/tan(theta*pi/180);
	double thetasp=2*gamasp+theta;
	double Ctheta=0;
	if(thetasp<7)
		Ctheta=(thetasp-7)/2;
	double D=0;
	D=-10*log10(1+2*sin(gamasp*pi/180)/(cos(thetasp*pi/180)*sin((gamasp+theta)*pi/180)));
	double pr=0;

	/*re1=re-hs;
	re2=re+ha;

	double fai=acos(((re+hs)*(re+hs)+re2*re2-dsa*dsa*1000000)/(2*(re+hs)*re2));//地心角,rad
	
	double xtr=fai*re;
	double dt1=hs*xtr/(ha+hs);
	double dr1=xtr-dt1;

	double fait=dt1/re;
	double fair=dr1/re;
	rsa1=sqrt(re*re+(re+hs)*(re+hs)-2*re*(re+hs)*cos(fait));
	rsa2=sqrt(re*re+(re+ha)*(re+ha)-2*re*(re+ha)*cos(fair));
	rsa=rsa1+rsa2;
	double dr=rsa-dsa*1000;

	double psi=(acos((rsa1*rsa1+re*re-pow(re+hs,2))/(2*re*rsa1))-pi/2)*1000;//擦地角，mrad*/

	CMultDlg Dlg;

	dxcomplex n2; 
	n2.real= epsilon;
	n2.imag=(60*sigma*wl);

	dxcomplex Rh;
	dxcomplex Rv;
	double rphase=0;
	double rmag=0;



	double cos2;
	double sin1;
	cos2 = cos(psi/1000) * cos(psi/1000);	
	sin1 = sin(psi/1000);
	double Rr = 0;
	double Ri = 0;

	dxcomplex temp1;
	dxcomplex temp2;


	n2.real=n2.real-cos2;
	Dlg.ontrt(&n2,2);	
	temp1.real=sin1-n2.real;
	temp2.real=sin1+n2.real;
	temp1.imag=-n2.imag;
	temp2.imag=n2.imag;

			
	Dlg.ocdiv(&temp1,&temp2);
	Rh = temp1;


	temp1.real=n2.real*sin1;
	temp1.imag=n2.imag*sin1;
	temp2.real=n2.real-cos2;
	temp2.imag=n2.imag-cos2;
	Dlg.ontrt(&temp2,2);
	Rv.imag=temp1.imag-temp2.imag;
	Rv.real=temp1.real-temp2.real;
	temp1.imag=temp1.imag+temp2.imag;
	temp1.real=temp1.real+temp2.real;
	Dlg.ocdiv(&Rv,&temp1);

	if (ipol==0)
	{
		Rr=Rh.real;
		Ri=Rh.imag;
		//if(f>0.3)
		//{
		//	Rr=-1;
		//	Ri=0;
		//}

	}
	else if(ipol==1)
	{
		Rr=Rv.real;
		Ri=Rv.imag;

	}
	else
	{
		Rr=(Rh.real+Rv.real)/2;
		Ri=(Rh.imag+Rv.imag)/2;

	}
	

	if (Rr > 0)
	{
		rphase = atan(Ri / Rr);
	} 
	else
	{
		if (Ri >= 0)
		{
			rphase = atan(Ri / Rr) + pi;
		}
		else
		{
			rphase = atan(Ri / Rr) - pi;
		}
	}
	rmag = sqrt(Rr * Rr + Ri * Ri);


	


	*dt=dr*10/3;
	*damp=20*log10(rmag);

	double pL=0;
	if(season==1)
		pL=20;
	if(season==2)
		pL=5;
	if(season==3)
		pL=5;
	if(season==4)
		pL=1;


	double k=pow(pL,1.5)*pow(10,-5.4);

	double epsp=fabs(ha-hs)/dsa;
	
	if(theta>5)
		*AMout=0;
	else
	    *AMout=10*log10(k*pow(dsa,3.3)*pow(f*1000,0.93)*pow(1+fabs(epsp),-1.1)*pow(psi,-1.2))-10*log10(100-pout);


	pr=20*log10(rmag)+Ctheta+D;
	double alpha=pow(10,pr/10)/(1+pow(10,pr/10));

	*AMout=5+10*log10(1+pow(10,pr/10));

	*Lr=-10*log10(1+rmag*rmag-2*rmag*cos(fko*dr+rphase));
}


void CMultDlg::CalAttIo(double f, double theta, double *attio)
{
	*attio=2500/(sin(theta*3.1415926/180)*f*f*pow(10,6));

}


void CMultDlg::CalFrady(double N[], double h[], double ht, double f, double theta, double fai, double *Frady, double *attfrady, double *XPD)
{
	double pi=3.1415926;
	theta=theta*pi/180;
	fai=fai*pi/180;
	ht=ht/1000;
	double z;

	double BL;
//	double lon,colat;
	double B1,B2,B3;
	static double Bx[1000]={29211.9,29141.3,29070.8,29000.6,28930.7,28860.9,28791.4,
		28722.0,28653.0,28584.1,28515.4,28447.0,28378.8,28310.8,28243.0,28175.4,
		28108.1,28040.9,27974.0,27907.3,27840.7,27774.5,27708.4,27642.5,27576.8,
		27511.3,27446.1,27381.0,27316.1,27251.5,27187.0,27122.8,27058.7,26994.9,
		26931.2,26867.8,26804.6,26741.5,26678.6,26616.0,26553.5,26491.2,26429.1,
		26367.3,26305.6,26244.0,26182.7,26121.6,26060.7,25999.9,25939.4,25879.0,
		25818.8,25758.8,25699.0,25639.3,25579.9,25520.6,25461.5,25402.6,25343.9,
		25285.3,25227.0,25168.8,25110.8,25052.9,24995.3,24937.8,24880.5,24823.4,
		24766.4,24709.6,24653.0,24596.6,24540.3,24484.2,24428.3,24372.5,24316.9,
		24261.5,24206.2,24151.1,24096.2,24041.4,23986.8,23932.4,23878.1,23824.0,
		23770.1,23716.3,23662.7,23609.2,23555.9,23502.8,23449.8,23397.0,23344.3,
		23291.8,23239.4,23187.2,23135.2,23083.3,23031.6,22980.0,22928.6,22877.3,
		22826.2,22775.2,22724.4,22673.7,22623.2,22572.8,22522.6,22472.5,22422.6,
		22372.8,22323.2,22273.7,22224.3,22175.1,22126.1,22077.1,22028.4,21979.8,
		21931.3,21882.9,21834.7,21786.7,21738.7,21691.0,21643.3,21595.8,21548.5,
		21501.2,21454.1,21407.2,21360.4,21313.7,21267.1,21220.7,21174.4,21128.3,
		21082.3,21036.4,20990.7,20945.0,20899.6,20854.2,20809.0,20763.9,20718.9,
		20674.1,20629.4,20584.8,20540.4,20496.0,20451.8,20407.8,20363.8,20320.0,
		20276.3,20232.7,20189.3,20146.0,20102.8,20059.7,20016.7,19973.9,19931.2,
		19888.6,19846.1,19803.8,19761.6,19719.4,19677.5,19635.6,19593.8,19552.2,
		19510.7,19469.3,19428.0,19386.8,19345.7,19304.8,19264.0,19223.3,19182.7,
		19142.2,19101.8};
	static double By[1000]={-898.7,-895.3,-891.9,-888.5,-885.2,-881.9,-878.6,-875.3,-872.0,
		-868.7,-865.5,-862.2,-859.0,-855.8,-852.6,-849.4,-846.3,-843.1,-840.0,-836.8,
		-833.7,-830.6,-827.5,-824.5,-821.4,-818.4,-815.3,-812.3,-809.3,-806.3,-803.3,
		-800.4,-797.4,-794.5,-791.5,-788.6,-785.7,-782.8,-780.0,-777.1,-774.2,-771.4,
		-768.6,-765.8,-763.0,-760.2,-757.4,-754.6,-751.9,-749.1,-746.4,-743.7,-741.0,
		-738.3,-735.6,-732.9,-730.2,-727.6,-724.9,-722.3,-719.7,-717.1,-714.5,-711.9,
		-709.3,-706.8,-704.2,-701.7,-699.1,-696.6,-694.1,-691.6,-689.1,-686.7,-684.2,
		-681.7,-679.3,-676.8,-674.4,-672.0,-669.6,-667.2,-664.8,-662.4,-660.1,-657.7,
		-655.4,-653.0,-650.7,-648.4,-646.1,-643.8,-641.5,-639.2,-637.0,-634.7,-632.5,
		-630.2,-628.0,-625.8,-623.5,-621.3,-619.1,-617.0,-614.8,-612.6,-610.5,-608.3,
		-606.2,-604.0,-601.9,-599.8,-597.7,-595.6,-593.5,-591.4,-589.4,-587.3,-585.2,
		-583.2,-581.2,-579.1,-577.1,-575.1,-573.1,-571.1,-569.1,-567.1,-565.1,-563.2,
		-561.2,-559.3,-557.3,-555.4,-553.5,-551.6,-549.7,-547.8,-545.9,-544.0,-542.1,
		-540.2,-538.4,-536.5,-534.7,-532.8,-531.0,-529.2,-527.3,-525.5,-523.7,-521.9,
		-520.1,-518.3,-516.6,-514.8,-513.0,-511.3,-509.5,-507.8,-506.1,-504.3,-502.6,
		-500.9,-499.2,-497.5,-495.8,-494.1,-492.4,-490.8,-489.1,-487.4,-485.8,-484.1,
		-482.5,-480.9,-479.2,-477.6,-476.0,-474.4,-472.8,-471.2,-469.6,-468.0,-466.4,
		-464.9,-463.3,-461.7,-460.2};
	static double Bz[1000]={39169.8,39071.6,38973.8,38876.3,38779.1,38682.2,38585.7,38489.5,
		38393.6,38298.0,38202.7,38107.8,38013.1,37918.8,37824.8,37731.1,37637.7,37544.6,
		37451.8,37359.3,37267.2,37175.3,37083.8,36992.5,36901.5,36810.9,36720.5,36630.4,
		36540.6,36451.1,36361.9,36273.0,36184.4,36096.1,36008.1,35920.3,35832.8,35745.6,
		35658.8,35572.1,35485.8,35399.7,35313.9,35228.4,35143.2,35058.2,34973.5,34889.2,
		34805.0,34721.1,34637.5,34554.2,34471.1,34388.3,34305.8,34223.5,34141.5,34059.8,
		33978.3,33897.1,33816.1,33735.4,33654.9,33574.7,33494.8,33415.1,33335.7,33256.5,
		33177.6,33098.9,33020.4,32942.2,32864.3,32786.6,32709.2,32631.9,32555.0,32478.3,
		32401.8,32325.5,32249.5,32173.8,32098.2,32022.9,31947.9,31873.0,31798.5,31724.1,
		31650.0,31576.1,31502.4,31429.0,31355.8,31282.8,31210.0,31137.5,31065.2,30993.1,
		30921.3,30849.6,30778.2,30707.0,30636.0,30565.3,30494.7,30424.4,30354.3,30284.4,
		30214.7,30145.3,30076.0,30007.0,29938.1,29869.5,29801.1,29732.9,29664.9,29597.1,
		29529.6,29462.2,29395.0,29328.1,29261.3,29194.7,29128.4,29062.2,28996.3,28930.5,
		28864.9,28799.6,28734.4,28669.5,28604.7,28540.1,28475.7,28411.5,28347.5,28283.7,
		28220.1,28156.7,28093.5,28030.4,27967.6,27904.9,27842.4,27780.1,27718.0,27656.1,
		27594.4,27532.8,27471.4,27410.2,27349.2,27288.4,27227.8,27167.3,27107.0,27046.9,
		26987.0,26927.2,26867.6,26808.2,26749.0,26689.9,26631.1,26572.4,26513.8,26455.4,
		26397.3,26339.2,26281.4,26223.7,26166.2,26108.8,26051.7,25994.6,25937.8,25881.1,
		25824.6,25768.3,25712.1,25656.0,25600.2,25544.5,25488.9,25433.6,25378.3,25323.3,
		25268.4};
	static float Hb[1000]={60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.0,115.0,
		120.0,125.0,130.0,135.0,140.0,145.0,150.0,155.0,160.0,165.0,170.0,175.0,180.0,
		185.0,190.0,195.0,200.0,205.0,210.0,215.0,220.0,225.0,230.0,235.0,240.0,245.0,
		250.0,255.0,260.0,265.0,270.0,275.0,280.0,285.0,290.0,295.0,300.0,305.0,310.0,
		315.0,320.0,325.0,330.0,335.0,340.0,345.0,350.0,355.0,360.0,365.0,370.0,375.0,
		380.0,385.0,390.0,395.0,400.0,405.0,410.0,415.0,420.0,425.0,430.0,435.0,440.0,
		445.0,450.0,455.0,460.0,465.0,470.0,475.0,480.0,485.0,490.0,495.0,500.0,505.0,
		510.0,515.0,520.0,525.0,530.0,535.0,540.0,545.0,550.0,555.0,560.0,565.0,570.0,
		575.0,580.0,585.0,590.0,595.0,600.0,605.0,610.0,615.0,620.0,625.0,630.0,635.0,
		640.0,645.0,650.0,655.0,660.0,665.0,670.0,675.0,680.0,685.0,690.0,695.0,700.0,
		705.0,710.0,715.0,720.0,725.0,730.0,735.0,740.0,745.0,750.0,755.0,760.0,765.0,
		770.0,775.0,780.0,785.0,790.0,795.0,800.0,805.0,810.0,815.0,820.0,825.0,830.0,
		835.0,840.0,845.0,850.0,855.0,860.0,865.0,870.0,875.0,880.0,885.0,890.0,895.0,
		900.0,905.0,910.0,915.0,920.0,925.0,930.0,935.0,940.0,945.0,950.0,955.0,960.0,
		965.0,970.0,975.0,980.0,985.0,990.0,995.0,1000.0};
//	lon=112.4;
//	colat=47.3;
	z=0;
		double tec=0;

	for (int i=0;h[i+1]<=ht;i++)
	{	
			for (int j=0;j<189;j++)
			{
				if (Hb[j]==h[i])
				{
					B1=Bx[j];
					B2=By[j];
					B3=Bz[j];
					break;
				}
			}
			BL=fabs(B1*cos(theta)*cos(fai)+B2*cos(theta)*sin(fai)-B3*sin(theta));
			z+=BL*N[i]*(h[i+1]-h[i]);	
			tec+=N[i]*(h[i+1]-h[i]);				

	}
			for (int j=0;j<188;j++)
			{
				if (ht>=h[j]&&ht<h[j+1])
				{
					B1=Bx[j];
					B2=By[j];
					B3=Bz[j];
					break;
				}
			}

			BL=fabs(B1*cos(theta)*cos(fai)+B2*cos(theta)*sin(fai)-B3*sin(theta));
			z+=BL*N[i]*(ht-h[i]);
			tec+=BL*N[i]*(ht-h[i]);

/*		hd=(ht-h[i])/2;
		Nd=(N[i+1]-N[i])*(ht-h[i])/((h[i+1]-h[i])*2);

	for (int j=0;j<2;j++)
	{
			hx=h[i]+j*hd;
			Nx=N[i]+j*Nd;
			BL=fabs(CMyDataFun::Fx(hx,lon,colat)*cos(theta)*cos(fai)
				+CMyDataFun::Fy(hx,lon,colat)*cos(theta)*sin(fai)
				-CMyDataFun::Fz(hx,lon,colat)*sin(theta));
			z+=BL*Nx*hd;		
	}
*/
//	z=z*1000*1E-9;
	z=z*1000*1E-9/sin(theta);
	z=(1.6E-19)*(1.6E-19)*(1.6E-19)*z/((3E8)*(8.85E-12)*(9.1E-31)*(9.1E-31)*4*3.1415926535*3.1415926535*f*f);


//		double test=(1.6E-19)*(1.6E-19)*(1.6E-19)/((3E8)*(8.85E-12)*(9.1E-31)*(9.1E-31)*4*3.1415926535*3.1415926535);

	*attfrady=-10*log10(cos(z)*cos(z));

	*XPD=-20*log10(tan(z));

	
	z=z*180/pi;
	*Frady=z;
		
}


void CMultDlg::CalSeSan(float N[], float h[], float ht, float f, float df, float theta, float *ShiYSS, float *XiangWSS, float *Pe)
{
	////////////////////////////////////////////////////////////////////////////////
	//
	//						N[1000]:电子浓度剖面
	//						h[1000]:与N[1000]相对应的高度;	
	//						ht:指定高度
	//						f:载频(hz)
	//						df:信号带宽(hz)
	//			            theta:地空水平仰角,度
    //
	//
	//					高度与折射率剖面的排列顺序为从60km处到高空
	//
	//		其中ht，h[1000]的单位为km,N[1000]单位为1/立方米	
	//			
	//					此函数用来到指定高度处的射线的时延色散ShiYSS与相位色散XiangWSS
	//
	//					
	//
	////////////////////////////////////////////////////////////////////////////////

//	float x;
	*ShiYSS=2.68E-7*CMultDlg::TEC(N,h,ht,theta)*df/(f*f*f);
	*XiangWSS=8.44E-7*df*CMultDlg::TEC(N,h,ht,theta)/(f*f);
    *Pe=pow(*ShiYSS,5.08)*(4.17E-7);
	//	return x;

}

void CMultDlg::CalShiY(float N[], float h[], float ht, float f, float theta, float *ShiY)
{
	////////////////////////////////////////////////////////////////////////////////
	//
	//						N[1000]:电子浓度剖面
	//						h[1000]:与N[1000]相对应的高度;	
	//						ht:指定高度
	//						f:载频(hz)
	//						theta:地空水平仰角,度
	//			
	//					高度与折射率剖面的排列顺序为从60km处到高空
	//
	//		其中ht，h[1000]的单位为km,N[1000]单位为1/立方米	
	//			
	//					此函数用来到指定高度处的射线的时延
	//
	//						计算结果为s
	//
	////////////////////////////////////////////////////////////////////////////////

//	float x;
	*ShiY=1.345E-7*CMultDlg::TEC(N,h,ht,theta)/(f*f);
	//	return x;

}


float CMultDlg::TEC(float N[], float h[], float ht, float theta)
{
	////////////////////////////////////////////////////////////////////////////////
	//
	//						N[1000]:电子浓度剖面
	//						h[1000]:与N[1000]相对应的高度;	
	//						ht:指定高度
	//			         	theta:地空水平仰角,度
	//
	//					高度与折射率剖面的排列顺序为从60km处到高空
	//
	//		其中ht，h[1000]的单位为km,N[1000]单位为1/立方米	
	//			
	//					此函数用来到指定高度处的电子总量
	//
	//						计算结果单位为		1/平方米
	//
	////////////////////////////////////////////////////////////////////////////////


	float x;
	x=0;
	float NN12;

theta=theta*3.1415926/180;
//	CSatdll sat;

	for (int i=0;i<999;i++)
	{
		if	(h[i]<ht&&h[i+1]<=ht)	x=x+(N[i]+N[i+1])*(h[i+1]-h[i])/2;
		if	(h[i]<=ht&&h[i+1]>ht)
		{
			NN12=N[i]+(ht-h[i])*(N[i+1]-N[i])/(h[i+1]-h[i]);
			x=x+(N[i]+NN12)*(ht-h[i])/2;
			break;
		}
	}

	
	x=x*1000;

	 //仰角的计算 调用动态连接库

	x=x/sin(theta);
	return x;
}

void CMultDlg::CalIoSin(double l, double fai, double fs, double hs, double f, double theta, double *Ais)
{
	////////////////////////////////////////////////////////////////////////////////
	//                  l=120*pi/180;              地面站经度109.52
	//                  fai=(0.0001:5:50)*pi/180;  地面站纬度18.23
	//                  fs=120*pi/180;    卫星定点位置176
	//                  hs=0.02;             地面站高度,km
	//                  f=2.45;          %频率，GHz

	//			            theta:地空水平仰角,度

	//			
	//					此函数用来计算电离层闪烁指数与衰减
	//
	//					
	//
	////////////////////////////////////////////////////////////////////////////////
	double pi=3.1415926;
	double a=4.8516;
	double lt=0;
	double fait=0;
	double faitm=0;
	double ltm=0;
	double lmp=-70.5*pi/180;
	double faimp=78.8*pi/180;
	double alpha=0;
	double faiq=0;
	double lq=0;
	double faiqm=0;
	double lqm=0;
	double beta=0;
	double w=0.71;
	double s4=0;
	double psi=0;

	l=l*pi/180;
	fai=fai*pi/180;
	fs=fs*pi/180;
	theta=theta*pi/180;




	lt=l-pi/2;
	psi=pi-atan(tan(fs-l)/sin(fai));     //卫星方位角

	faitm=asin(sin(fait)*sin(faimp)+cos(fait)*cos(faimp)*cos(lt-lmp));
	ltm=pi-asin(cos(fait)*sin(lt-lmp)/cos(faitm));

	alpha=pi/2-theta-asin((6370+hs)*cos(theta)/6670);
	faiq=asin(sin(fai)*cos(alpha)+cos(fai)*sin(alpha)*cos(psi));
	lq=l+asin(sin(alpha)*sin(psi)/cos(faiq));
	faiqm=asin(sin(faiq)*sin(faimp)+cos(faiq)*cos(faimp)*cos(lq-lmp));
	lqm=pi-asin(cos(faiq)*sin(lq-lmp)/cos(faiqm));

	beta=abs(lqm-ltm);

	w=0.71;

	s4=a/(pow(f,1.5)*sqrt(sin(theta))*exp(beta/w));
	//*Ais=7.1*s4*s4+12.59*s4-0.25;
	*Ais=27.5*pow(s4,1.26)/sqrt(2);

}


int CMultDlg::GetFileData(CString sPath, double** fData1, double** fData2)
{
	// 打开文件
	CFile file;
	file.Open(sPath, CFile::modeReadWrite);
	
	// 得到文件长度
	DWORD dwLen = file.GetLength();
	
	// 根据文件长度分配内存
	char* pData = new char [dwLen];
	memset(pData, 0, sizeof(char) * dwLen);
	
	// 将文件里的数据读到char型数组
	file.Read(pData, dwLen);
	
	// 关闭文件
	file.Close();
	
	// 统计行数
	int nRow = 1;
	for (DWORD i = 0; i < dwLen; i++)
	{
		if (pData[i] == 13)
			nRow++;
	}
	
	// 初始化double数组
	*fData1 = new double[nRow];
	memset(*fData1, 0, sizeof(double) * nRow);
	
	// 初始化float2数组
	*fData2 = new double[nRow];
	memset(*fData2, 0, sizeof(double) * nRow);
	
	// 初始化临时数组
	char temp[256];
	memset(temp, 0, sizeof(temp));
	
	// temp数组索引变量
	int j = 0;
	
	// fData1数组索引变量
	int k = 0;
	
	// fData2数组索引变量
	int l = 0;
	
	// 数据分段提取
	for (i = 0; i < dwLen; i++)
	{
		// 检测到回车符后将当前保存到temp的字符转换为double变量
		if (pData[i] == 10)
			continue;
		else if (pData[i] == ' ' || pData[i] == '\t' )	
		{
			*(*fData1 + k) = (double)atof(temp);
			
			// 清空temp数组以备下一个数使用
			memset(temp, 0, sizeof(temp));
			j = 0;
			k++;
		}
		else if (pData[i] == 13)
		{
			*(*fData2 + l) = (double)atof(temp);
			
			// 清空temp数组以备下一个数使用
			memset(temp, 0, sizeof(temp));
			j = 0;
			l++;
		}		
		// 检测到非回车符和换行符的内容当作数据处理，保存到temp数组
		else
		{
			temp[j] = pData[i];
			j++;
		}
	}
	
	if (j > 0)
		*(*fData2 + l) = (double)atof(temp);
	else
		nRow--;
	
	// 释放内存
	delete[] pData;
	
	return nRow;
}