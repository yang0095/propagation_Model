// TropoScatter.cpp: implementation of the CTropoScatter class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "TropoScatter.h"
#include <math.h>


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTropoScatter::CTropoScatter()
{

}

CTropoScatter::~CTropoScatter()
{

}



double CTropoScatter::GetAe(double dN)
{
	//得到等效地球半径 ae(km)

	double A=6371.0;//实际地球半径(km)
	double k=157/(157-dN);
	double ae=k*A;
	
	return ae;

}
/*
参数说明：
freq:频率, MHz 
hts,hrs: 发射、接收天线高度，m
d: 距离，km
p:时间百分比，%
Gt,Gr:天线增益，dB
dN：近地面1km折射率梯度
climate:气候带, 气候带：1/2/3/4/5/6/'7a'->7, '7b'->8
los:损耗(dB)
ipath：true-超视距路径，false-视距路径

*/

void CTropoScatter::TropoScatCal(double freq, double hts, double hrs, double d, double p,  double Gt, double Gr,double dN, int climate, double &los, bool &ipath)
{
	double ae=GetAe( abs(dN) );

	//计算散射角的程序
	double Theta0, tht, thr;
	tht = -1*acos(ae/(ae+hts*1.e-3)) *1e3;
	double thtd=(hrs-hts)/d-1e3*d/2/ae;

	if (tht>thtd) //超视距路径
	{
		ipath=true;

		thr = -1*acos(ae/(ae+hrs*1.e-3)) *1e3;
		double the=d*1e3/ae;
		Theta0=the+tht+thr;
		double sigma;
		LosTroposcatter(freq,d,p,Theta0,Gt,Gr,ae,climate,los,sigma); //(zr)Theta0重新计算
	} 
	else //视距路径
	{
		ipath=false;
		los = 0;
	}

}

void CTropoScatter::LosTroposcatter(double freq, double d, double p,  double Theta0, double Gt, double Gr,double ae, 
									int climate,double &los, double &sigma)
{
	//计算对流层散射损耗，用张院士的公式
	//把方差sigma也一起输出
	
	//参数说明：freq:频率(Mhz), Gt,Gr:天线增益，p:时间百分比，climate:气候带
	//    输出：los:损耗(dB), sigma:方差
	
	//by zhangr		
	
	//step 2
	double M,gama;
	GetPara(climate,M,gama);
	
	//step 4
	double H=1e-3*Theta0*d/4;
	double h=1e-6*Theta0*Theta0*ae/8;
	double LN=20*log10(5+gama*H)+4.34*gama*h;
	
	//step 5
	double Cq=-Qi(p/100)/1.2816;
	double Y_90=Y90(freq,h,climate,ae,Theta0);
	double Yq=Cq*Y_90;
	
	//step 6
	double Lc=0.07*exp(0.055*(Gt+Gr)); //天线口面介质耦合损耗
	
	los=M+30*log10(freq)+10*log10(d)+30*log10(Theta0)+LN+Lc-Gt-Gr-Yq;
	
	
	sigma=-1*Y_90/1.2816;

}

void CTropoScatter::GetPara(int climate, double &M, double &gamma)
{
	//根据气候带获取气象参数
	switch (climate)
	{
	case 1:
		M=39.60;
		gamma=0.33;
		break;
		
	case 2:
		M=29.73;
		gamma=0.27;
		break;
		
	case 3:
		M=19.30;
		gamma=0.32;
		break;
		
	case 4:
		
	case 5:
		M=38.50;
		gamma=0.27;
		break;
		
	case 6:
		M=29.73;
		gamma=0.27;
		break;
		
	case 7: //7a
		M=33.20;
		gamma=0.27;
		break;
		
	case 8: //7b
		M=26.00;
		gamma=0.27;
		break;
		
	default:
		break;
		
	}

}

double CTropoScatter::Qi(double x)
{
	//累积正态分布的逆函数
	x=max(1.e-6,x);
	x=min(0.999999,x);
	
	double y;
	if (x<=0.5)
	{
		y=T(x)-Xi(x);
	} 
	else
	{
		y=Xi(1-x)-T(1-x);
	}
	
	return y;

}



double CTropoScatter::Y90(double freq, double h, int climate, double ae, double Theta0)
{
	//计算Y(90)：freq(MHz), h(km), climate(气候带)
	//7a:7,  7b:8
	
	double y;
	
	double ds=ae*Theta0/1000;
	
	//气候带1
	const int n1=9;
	double d1[n1]={100,200,300,400,500,600,700,800,900};
	double y1[n1]={-8.1,-6.6,-5.2,-4.5,-4,-4,-3.4,-3.2,-3.1};
	
	//气候带3
	const int n3=11;
	double d3[n3]={100,150,200,250,300,360,400,450,500,550,900};
	double y3[n3]={-11.0,-12.3,-13.0,-12.5,-11.5,-10.0,-9.2,-8.8,-8.6,-8.5,-8.5};
	
	//气候带4
	const int n4=7;
	double d4[n4]={100,200,300,400,500,600,900};
	double y4[n4]={-11.5,-9.8,-7.6,-5.9,-4.3,-4,-4};	
	
	switch (climate)
	{
	case 2:
	case 6:
	case 7:
		y=-2.2  -(8.1-2.3e-4*freq) *exp(-0.137*h);
		break;
	case 8:
		y=-9.5  -3*exp(-0.137*h);
		break;
	case 1:
		Interp1(d1,y1,n1,&ds,&y,1,"nearest");
		break;		
	case 3:
		Interp1(d3,y3,n3,&ds,&y,1,"nearest");
		break;
	case 4:
		Interp1(d4,y4,n4,&ds,&y,1,"nearest");		
		break;	
		
	}
	
	return y;

}

double CTropoScatter::T(double x)
{
	//累积正态分布的逆函数__T(x)
	
	double y=sqrt(-2*log(x));
	return y;
}

double CTropoScatter::Xi(double x)
{
	//累积正态分布的逆函数__Xi(x)
	
	double c0=2.51551698;
	double c1=0.802853;
	double c2=0.010328;
	double d1=1.432788;
	double d2=0.189269;
	double d3=0.001308;
	
	double fenzi=( (c2*T(x)+c1)*T(x) )  +c0;
	double fenmu=( (d3*T(x)+d2)*T(x) +d1 ) *T(x) +1;
	
	double y=fenzi/fenmu;
	
	return y;
}



void CTropoScatter::Interp1(double *x, double *y, int n, double *x1, double *y1, int n1, CString method)
{
	//对数组的插值运算，参照matlab的方法
	
	//输入:x,y:待插值函数；x1,y1:输出的插值函数，n:x1,y1的大小	
	//method: linear表示数据超界用线性插值，nearest表示数据超界用两端数
	
	double h0,h1,h2; 
	double att0,att1,att2;
	
	
	int j=1;
	int k=0;
	for (int i=0;i<n1;i++) 
	{
		h0=x1[i];
		
		if (h0<=x[0]) //小于最小
		{
			if (method=="linear")
			{
				y1[i]=Linear( x[0],y[0],x[1],y[1],x1[i] );
			}
			else if (method=="nearest")
			{
				y1[i]=y[0];
			}			
		} 
		else if (h0>=x[n-1]) //大于最大
		{
			if (method=="linear")
			{
				y1[i]=Linear( x[n-2],y[n-2],x[n-1],y[n-1],x1[i] );
			}
			else if (method=="nearest")
			{
				y1[i]=y[n-1];
			}
			
		}
		else //中间
		{
			
			h1=x[k];
			h2=x[j];
			
			while (h0>h2)
			{
				j=j+1;
				k=j-1;
				
				h1=x[k];
				h2=x[j];
			}
			
			att1=y[k];
			att2=y[j];
			att0=att1+(att2-att1)*(h0-h1)/(h2-h1);
			
			y1[i]=att0;
			
		}		
	}	
}

double CTropoScatter::Linear(double x1, double y1, double x2, double y2, double x0)
{
	//一维线性插值，可以进行重载
	double k=(y2-y1)/(x2-x1);
	double y0=k*(x0-x1)+y1;
	return y0;
	
}





















