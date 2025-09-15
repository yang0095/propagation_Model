// P1812_Fun.cpp: implementation of the CP1812_Fun class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "P1812.h"
#include "P1812_Fun.h"
#include <math.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CP1812_Fun::CP1812_Fun()
{

}

CP1812_Fun::~CP1812_Fun()
{

}


/*
����˵����
P1812-6������������ITU��վ�ϵ�һ�£�(zr)2023-4-4: Ӧ�û�����ЩС��ͬ�ģ�

����˵����
freq����Ƶ�ʣ�GHz
p����ʱ��ٷֱȣ�p%
pL�����ص�ٷֱȣ�pL%
latMid����·���е��γ�ȣ���
htg��������վ�����߸߶ȣ�agl. m
hrg��������վ�����߸߶ȣ�agl. m
pol����������0:ˮƽ,1:��ֱ
erp������Ч���书�ʣ���kW����λ�ȽϺ��ʣ�(kW)

n����n>=3
dn����dn�ļ����һ���ǵȼ��ģ�ע�ⵥλ��km
hn����m
Rn�������ӵ���߶ȣ�m
climn�����������Ϣ�������P452��һ�£�Ҳ�������ԭ�������еı�����Coastal land A1(1), Inland A2 (2) or Sea B (3);
dN����������1km���������ݶȣ�dN>0����Ҫȡ��ֵ��
N0����

Lb�����������	
Ep���� ���ճ�ǿ
errorMes����������Ϣ
*/
bool CP1812_Fun::P1812_6_ITUwb(double freq, double p, double pL, double latMid, double htg, double hrg, int pol, double erp,//�������е�����
							   int n, double *dn, double *hn,double *Rn, int *climn, double dN, double N0, //�߳����ݺͻ�������
							   double &Lb,double &Ep, CString &errorMes) //�������					  
							 
{
	//����������(zr)?����������ʱ���ټӣ�	
	if (freq<0.03 || freq>6)
	{
		errorMes="Ƶ�ʷ�ΧΪ[0.03, 6]GHz��";
		return false;
	}
	if (p<1 || p>50)
	{
		errorMes="ʱ����ʷ�ΧΪ[1, 50]%��";
		return false;
	}

	if (pL<1 || pL>99)
	{
		errorMes="�ص���ʷ�ΧΪ[1, 99]%��";
		return false;
	}

	if (latMid<-80 || latMid>80)
	{
		errorMes="����·���е�γ�ȷ�ΧΪ[-80, 80]�ȣ�";
		return false;
	}

	if (htg<=0 || htg>3000)
	{
		errorMes="�������ߵĸ߶ȷ�ΧΪ(0, 3000]m��";
		return false;
	}

	if (hrg<=0 || hrg>3000)
	{
		errorMes="�������ߵĸ߶ȷ�ΧΪ(0, 3000]m��";
		return false;
	}

	if (pol!=0 && pol!=1)
	{
		errorMes="���߼�����ȡֵΪ0��1��";
		return false;
	}
	if (n<3)
	{
		errorMes="·�������������Ϊ3����";
		return false;
	}
	
	//����Ԥ����
	double *gn=new double[n];
	*gn = *hn;
	for (int i=1;i<n-1;i++) //����˺ͽ��ն˲��Ӹ��ӵ���߶�
	{
		*(gn+i)=*(hn+i) + *(Rn+i);
	}
	*(gn+n-1) =*(hn+n-1);

	
	//��������
	dN = fabs(dN);
	double ae,b0,tau,w; 
	ParaInit(latMid,n,dn,climn,dN,ae,b0,tau,w);

	//ʱ��仯
	double Lb0p,Lbc;
	CompLossTime(freq,p,pol,htg,hrg,n,dn,hn,gn,climn,N0,ae,b0,tau,w,Lb0p,Lbc);

	//�ص�仯��
	Lb = CompLossLoc(freq,pL,hrg,*(Rn+n-1),Lb0p,Lbc);
		
		
	//���ת��Ϊ��ǿ
	Ep=L2E(freq,Lb);	
	Ep=Ep+10*log10(erp);


	//ɾ������
	delete []gn;

	return true;
}


/*
�������ܣ������ʼ������

����˵����
latt/lont/latr/lonr������/����վ�ľ�γ�ȣ�rad

���������ae,b0,tau,w;
*/
void CP1812_Fun::ParaInit(double latMid,int n,double *dn, int *climn, double dN,
							 double &ae,double &b0,double &tau,double &w)
{
	
	//��ȡdtm/dlm
	double dtm,dlm;
	GetDtmDlm(n,dn,climn,dtm,dlm);

	//(zr)����·�����ȱ����������á�w=1-dtm/d���������ǲ��Եģ��Ѹģ�
	w=GetW(n,dn,climn);	


	//��Ч����뾶ae
	ae=GetAe(dN);	
	
	//��ȡb0,tau
	GetB0(dlm,dtm,latMid,b0,tau);//ab0=A*3.0;	

}


/*
�������ܣ���dtm/dlm����Ҫ�ȼ����������Ϣ

���������
climn�����������Ϣ�������P452��һ�£�Ҳ�������ԭ�������еı�����Coastal land A1(1), Inland A2 (2) or Sea B (3);


���������
dtm:longest continuous land(inland+coastal/A2+A1),km
dlm:longest continuous inland(A2),km
  
˵��������Ժʿ����Ĳ�������������Ұ������ûȡ��dlm��ȫ·����ȡ����������dtm��ȡ
*/
void CP1812_Fun::GetDtmDlm(int n,double *dn,int *climn,double &dtm, double &dlm)
{
	
	
	int i1,i2;//����
	double dtmi,dlmi;
	
	dtm=0;
	i1=-1; i2=-1;
	for (int i=0;i<n;i++)
	{
		if (climn[i]!=3 && i1==-1 && i2==-1)
		{
			i1=i;
		}
		
		if ( climn[i]!=3 && i2==-1 && (i==n-1 || climn[i+1]==3) )
		{
			i2=i;
		}
		
		if (i1!=-1 && i2!=-1)
		{
			dtmi=dn[i2]-dn[i1];
			dtm=max(dtm,dtmi);
			
			i1=-1;i2=-1;
		}
	}
	
	dlm=0;
	i1=-1;i2=-1;
	for (i=0;i<n;i++)
	{
		if (climn[i]==2 && i1==-1 && i2==-1)
		{
			i1=i;
		}
		
		if (climn[i]==2 && i2==-1 && (i==n-1 || climn[i+1]!=2))
		{
			i2=i;
		}
		
		if (i1!=-1 && i2!=-1)
		{
			dlmi=dn[i2]-dn[i1];
			dlm=max(dlm,dlmi);
			
			i1=-1;i2=-1;
		}
	}	
	
}


/*
�������ܣ����ݻ������ͣ��õ�����·�����ȱ���
�õ�������ʽ���±�����

climn�����������Ϣ�������P452��һ�£�Ҳ�������ԭ�������еı�����Coastal land A1(1), Inland A2 (2) or Sea B (3);
	
(zr)����
*/
double CP1812_Fun::GetW(int n, double *dn, int *climn)
{

	double lw=0;
	if (climn[0]==3)
	{
		lw+=0.5*( dn[1]-dn[0] );
	}

	for (int i=1;i<n-1;i++)
	{
		if (climn[i]==3)
		{
			lw+=0.5*( dn[i+1]-dn[i-1] );
		}
	}

	if (climn[n-1]==3)
	{
		lw+=0.5*( dn[n-1]-dn[n-2] );
	}

	double w=lw/dn[n-1];

	return w;

}



/*
����˵�����õ���Ч����뾶 ae(km)
dNΪ��ֵ
*/
double CP1812_Fun::GetAe(double dN)
{
	
	double A=6371.0;//��Ч����뾶(km)

	double k=157/(157-dN);
	double ae=k*A;
	
	return ae;
	
}

/*
�������ܣ�
�õ�b0��ֵ����p!=50��ʱ�����Ч��������Ĳ�����Ҫ���ߵ�������Ĳ���

  
���������
dlm(km) dtm(km):��Ҫ�����ݵõ���
latMid����·���е�γ��(��)��
*/
void CP1812_Fun::GetB0(double dlm, double dtm, double latMid,double &b0,double &tau)
{
	
	tau=1-exp(-4.12e-4*pow(dlm,2.41));
	double umid=pow(10,-1*dtm/(16-6.6*tau))+pow(10,-5*(0.496+0.354*tau));
	double u1=pow(umid,0.2);
	
	double phi=latMid;//����ת��Ϊ��
	
	double u4;
	if (phi<=70)
	{
		u4=pow(u1,-0.935+0.0176*phi);
		b0=pow(10,-0.015*phi+1.67)*u1*u4;
	} 
	else
	{
		u4=pow(u1,0.3);
		b0=4.17*u1*u4;
	}
	
}


/*
����˵������Ե�ļ���������main��������ĺͳ�ǿ

����˵����
freq����Ƶ��(GHz)
p����ʱ��ٷֱ�
pL����(zr)���������
pol����������0:ˮƽ,1:��ֱ
htg/hrg��������ͽ������߸߶ȣ�agl.
n/dn/hn/gn/ngc����������km, �߶���m!
erp������Ч���书��


N0�����ر�������
ae������Ч����뾶
b0�����м����
tau�����м����
w��������·������

lb0p��·�����
lbc����ǿ
*/
void CP1812_Fun::CompLossTime(double freq,double p, int pol,double htg,double hrg,int n,double *dn, double *hn,double *gn, int *climn,
						 double N0,double ae,double b0,double tau,double w,
						 double &lb0p,double &lbc)
{
	double hts=*hn+htg;
	double hrs=*(hn+n-1)+hrg;
	double htc=hts;	    //(zr)�Ķ���//hts/hrs=htc/hrc 
	double hrc=hrs;
	double d=*(dn+n-1);
	
	//·���������	
	bool bpath;
	int ilt,ilr;
	double thetat,thetar,theta; 
	PathAnalyze(freq,hts,hrs,ae,n,dn,hn, bpath,ilt,ilr,thetat,thetar,theta);


	//________���ɿռ䴫�����,����ʱ�����������_________//
	double lbfs,lb0b;//lb0p	
	double dlt,dlr;
	dlt=*(dn+ilt);dlr=d-*(dn+ilr);
	lbfs=LosFreespace(freq,d,50,dlt,dlr);    //50%
	lb0p=LosFreespace(freq,d,p,dlt,dlr);     //p%
	lb0b=LosFreespace(freq,d,b0,dlt,dlr);    //b0%	


   
	//________�������+���ɿռ䴫�����_________//
	double ldp,ld50;
	//Ld(freq,htc,hrc,pol,p,b0,ae,w,n,dn,gn,ld50,ldp);	
	Ld1812_3(freq,htc,hrc,pol,p,b0,ae,w,n,dn,hn,gn,ld50,ldp);
	double lbd50,lbd;
	lbd50=lbfs+ld50;
	lbd=lb0p+ldp;


	//________������ɢ�䴫�����_________//
	//ʱ�����0.001%~50%��ͨ�õ�
	double lbs=LosTropoScatter(freq,p,d, theta, N0);
	
	
	//________�����������_________//
	double lba;		
	lba = LosDuct(freq,p,b0,hts,hrs,n,dn,hn,climn,ae,w,tau,ilt,ilr,thetat,thetar); 
	
	
	 //(zr)tothis�� 
	//________����p%ʱ���50%�ص�����,�������Ӳ���ģ���4.6_________//
	double Fi=GetFi(p,b0);
	double Lminb0p;
	if (p<b0)
	{
		Lminb0p=lb0p+(1-w)*ldp;
	} 
	else
	{
		Lminb0p=lbd50+(lb0b+(1-w)*ldp-lbd50)*Fi;//(zr)������ԭ��û��*ldp
	}

	double Lminbap;
	double yita=2.5;
	Lminbap=yita*log( exp(lba/yita)+exp(lb0p/yita) );

	double Lbda;
	double Fk=GetFjk(20,0.5,d);
	if (Lminbap>lbd)
	{
		Lbda=lbd;
	} 
	else
	{
		Lbda=Lminbap+(lbd-Lminbap)*Fk;
	}

	double Lbam;
	double Fj=GetFjk(0.3,0.8,theta);
	Lbam=Lbda+(Lminb0p-Lbda)*Fj;

	lbc=-5*log10( pow(10,-0.2*lbs)+pow(10,-0.2*Lbam) );

}

/*
�������ܣ�
��ȡ·���������������ʱ������ֻ��hn

˵����(zr)��ʱ��ĸģ���ʽ��ITUһ�±ȽϺ���???
*/
void CP1812_Fun::PathAnalyze(double freq,double hts,double hrs,double ae,int n,double *dn, double *hn,
								bool &bpath, int &ilt, int &ilr, double &thetat, double &thetar, double &theta)
{

	double di,hi,thi; 
	double thimax,thjmax,thtd;//�ֱ�Ϊ����ͽ��մ�����󴫲���
	double mid0=1.e3/2/ae;//�м����
	double d=*(dn+n-1);


	thtd=(hrs-hts)/d-mid0*d;

	//��ֵ��һ����С������-999.9
	thimax=-999;
	for (int i=1;i<n-1;i++)
	{
		di=*(dn+i);
		hi=*(hn+i);
		thi=(hi-hts)/di-mid0*di;

		if (thi>thimax)
		{
			thimax=thi; 
			ilt=i;
		}
	}


	if (thimax>thtd) //�ж��Ƿ�Ϊ���Ӿ���·  //���Ӿ�·��
	{
		bpath=true;

		thjmax=-999;
		for (i=ilt;i<n-1;i++) //������Ժʿ�Ĳ�һ����ilt~pPara.n-1
		{
			di=*(dn+i);
			hi=*(hn+i);
			thi=(hi-hrs)/(d-di)-mid0*(d-di);  
						
			if (thi>thjmax)
			{
				thjmax=thi;
				ilr=i;
			}
		}	
		
		thetat=thimax;		
		thetar=thjmax;

	}
	else  //�Ӿ�·��
	{
		bpath=false;
		
		ilt=Im50(freq,hts,hrs,n,dn,hn,ae);
		ilr=ilt;
		
		thetat=thtd;
		thetar=(hts-hrs)/d-mid0*d;
	}

	theta=thetat+thetar+mid0*2*d; //�Ǿ���


	/*
	//��ֵ���ṹ�����
	thtr.bPath=bpath;

	thtr.ilt=ilt;
	thtr.ilr=ilr;

	thtr.thetat=thetat;
	thtr.thetar=thetar;
	thtr.theta=theta;	
	*/

}


/*
����˵��������·���ϵ���ߵ�	
freq:GHz; �߶�:m; ����:km
*/
int CP1812_Fun::Im50(double freq,double h1,double h2, int n, double *dn, double *hn, double ae)
{
	
	
	double d1,d2,hi,Hi;
	double zeta,nu,numax;
	double lamda=0.3/freq;
	double d=*(dn+n-1);
	int im;
	
	
	//�������濼��0-n
	numax=-999; //��ʼֵ������С����
	for (int i=1;i<n-1;i++)
	{
		d1=*(dn+i);
		d2=d-d1;
		hi=*(hn+i);
		Hi=hi +1.e3*d1*d2/2/ae-(h1*d2+h2*d1)/d;
		zeta=cos(atan( 1.e-3*(h2-h1)/d ));
		nu=zeta*Hi*sqrt( 2.e-3*d/(lamda*d1*d2) );
		
		if (nu>numax)//�ж�
		{
			im=i;
			numax=nu;
		}		
	}
	
	return im;
}

/*
�������ܣ�
�������ɿռ����,�����������

���������freq:Ƶ��(GHz),d:·������(km),p:ʱ����ʣ�dlt,dlr:(km)
*/
double CP1812_Fun::LosFreespace(double freq,double d,double p,double dlt,double dlr)
{
	double lbfs=92.44+20*log10(freq)+20*log10(d);
	double esp=2.6*(1-exp(-(dlt+dlr)/10.0))*log10(p/50.);
	double lb0=lbfs+esp;
	
	return lb0;
}


/*
�������ܣ�����1812-3Э������������,��-6�汾û�иĶ���
(zr)����Ϊ����ʤ��д
���ld50��ldp

(zr)dn,gn��·������������������õĲ�ͬ��
*/

void CP1812_Fun::Ld1812_3(double freq,double htc,double hrc,int pol,double p,double b0,double ae,double w,int n,double *dn, double *hn, double *gn,
							 double &ld50, double &ldp)
{	
	double A=6371.0;//��Ч����뾶(km)

	double ap=ae;	
	double ld=LdDeltaBullington( freq, htc, hrc, pol, ap, w, n, dn, hn, gn);
	
	ld50=ld;
	if(p==50)
	{
		ldp=ld50;
	}
	else
	{ 
		ap=3.0*A;		
		ld=LdDeltaBullington( freq, htc, hrc, pol, ap, w, n, dn, hn, gn);

		double ldb=ld;
		double Fi=GetFi(p,b0);
		ldp=ld50+(ldb-ld50)*Fi;	
	}
	
}

/*
�������ܣ�
DeltaBullington����������ģ�ͼ���
*/
double CP1812_Fun::LdDeltaBullington(double freq,double htc,double hrc,int pol,double ap,double w,
										int n,double *dn,double *hn,double *gn)
{
	//������ĵ���Ҫ�������
	
	double d=*(dn+n-1);
	double lamda=0.3/freq;


	double hst,hsr; //ע�⣺�����õ���hn
	DeriSmoothEarthsurface(n,dn,hn, hst,hsr);
	double hstd,hsrd;
	DeriDiffPar(hst,hsr,htc,hrc, n,dn,hn, hstd, hsrd);
	double htesph=htc-hstd;
	double hresph=hrc-hsrd;
	
	
	//Bullington��һ�μ���
	double lbulla=LdBullington(freq,htc,hrc,n,dn,gn,ap);  //ע�⣺�����õ���gn

	
	//Bullington�ڶ��μ��㣬���߸߶ȱ��޸ģ�������߶���0
	double *gn0=new double[n];
	for (int i=0;i<n;i++)//�൱��hgc=0
	{
		*(gn0+i)=0;
	}
	double lbulls=LdBullington(freq,htesph,hresph,n,dn,gn0,ap);  //ע�⣺�����õ���gn0
	
	//���ε���������ļ���
	double ldsph=LdSphericalEarth( freq, pol, htesph, hresph, d, ap, w);
	
	
	//��������delta-Bullington���̣��õ�Ld
	double ld;
	double ld0=ldsph-lbulls;
	if(ld0>0)
	{
		ld=lbulla+ld0;
	}
	else
	{
		ld=lbulla;
	}//delta-Bullington���
	
	
	delete []gn0;
	
	return ld;
}

/*
����˵����bullington������ģ������ڲ������������Ԥ��

����˵����
  
���ߺ����ڣ�liuys
	
����˵�������ܲο���zhangr�ĳ���
*/
double CP1812_Fun::LdBullington(double freq,double htc,double hrc,int n,double *dn,double *gn,double ae)
{
	double luc;
	double d=*(dn+n-1);
	double lamda=0.3/freq;
	double stim,srim,v,dbp,lbull;	
	stim=Stim(htc,n,dn,gn,ae);
	double str=(hrc-htc)/d;

	if(stim<str)//LOS·��
	{
		v = Vmax(freq,htc,hrc,n,dn,gn,ae);//��ӵĺ���
		luc = LosKnifeEdge(v);
	}
	else 
	{
		srim = Srim(hrc,n,dn,gn,ae);
		dbp=(hrc-htc+srim*d)/(stim+srim);
		v = (htc+stim*dbp-(htc*(d-dbp)+hrc*dbp)/d)*sqrt(0.002*d/lamda/dbp/(d-dbp));
		luc = LosKnifeEdge(v);
	}

	lbull=luc+(1-exp(-luc/6))*(10+0.02*d);

	return lbull;
}

/*
����˵�����⻬����������� ��any distance������LOS��NLOS������4.3.2  ���͵����������

����˵����
  
���ߺ����ڣ�
	
����˵����
*/
double CP1812_Fun::LdSphericalEarth(double freq,int pol,double htesph,double hresph, double d,double ae,double w)
{
	double PI=3.14159265358979;//����pi

	double ldsph;
	double adft,ldft;
	double lamda=0.3/freq;

	double ap=ae;//apȡֵ������4.3.5
	double dlos=sqrt(2*ap)*(sqrt(0.001*htesph)+sqrt(0.001*hresph));
	if(d<dlos)
	{
		double mc=250*d*d/ap/(htesph+hresph);
		double c=(htesph-hresph)/(htesph+hresph);
		double b=2*sqrt((mc+1)/3/mc)*cos(PI/3+acos(3*c*sqrt(3*mc/pow(mc+1,3))/2)/3);
		double dse1=d*(1+b)/2;
		double dse2=d-dse1;

		double hse=(htesph-500*dse1*dse1/ap)*dse2/d+(hresph-500*dse2*dse2/ap)*dse1/d;
		double hreq=17.456*sqrt(dse1*dse2*lamda/d);

		if(hse>hreq)
		{ldsph=0;}
		else 
		{  
			double aem=500*pow(d/(sqrt(htesph)+sqrt(hresph)),2);
			adft=aem; 
			ldft=LdFirstTermComp(freq,d,htesph,hresph,adft,w,pol);

			if(ldft<0)
			{
				ldsph=0;
			}
			else
			{
				ldsph=(1-hse/hreq)*ldft;
			}
		}
	}
	else //����4.3.3������adft=ap����Ldft��������Ldsph����Ldft
	{
		adft=ap;
		ldft=LdFirstTermComp(freq,d,htesph,hresph,adft,w,pol);
		ldsph=ldft;
	}
	return ldsph;
}


//����������õĲ�ֵ���ӣ�����ʽ(30)
double CP1812_Fun::GetFi(double p, double b0)
{	
	double Fi;
	if(p==50.)
	{
		Fi=0;
	}
	else if (p>b0)
	{
		Fi=I(p/100)/I(b0/100);
	} 
	else
	{
		Fi=1;
	}
	
	return Fi;
}


/*
�������ܣ��Ƶ��⻬����������,�õ���һ���⻬�������(deriving the smooth-Earth surface);
          �ο�P1812-3; ��P1812-1��ͬ;
by zhangr, 2023-1-13
*/
void CP1812_Fun::DeriSmoothEarthsurface(int n, double *dn, double *hn, double &hst, double &hsr)
{	
	double d=*(dn+n-1);
	double v1=0;
	double v2=0;

	for (int i=1;i<n;i++)
	{
		v1 += (dn[i]-dn[i-1])*(hn[i]+hn[i-1]);
		v2 += (dn[i]-dn[i-1])*( hn[i]*(2*dn[i]+dn[i-1]) + hn[i-1]*(dn[i]+2*dn[i-1]) );	
	}
	hst = (2*v1*d-v2)/d/d;
	hsr = (v2-v1*d)/d/d;

}



/*
�������ܣ�
·����������� ����ģ�Ͳ�������ITU�Ľ��һ�£�
*/
void CP1812_Fun::DeriDuctPar(int n, double *dn,double *hn, double hts, double hrs, int ilt, int ilr, double hst, double hsr,
								  double &hm, double &hte,double &hre)
{
	hst=min(hst,*hn);
	hsr=min(hsr,*(hn+n-1));
	
	double d=*(dn+n-1);
	double m=(hsr-hst)/d;
	
	//hm
	double htemp;
	hm=-999;
	for (int i=ilt;i<=ilr;i++)
	{
		htemp=hn[i]-(hst+m*dn[i]);
		if (htemp>hm)
		{
			hm=htemp;
		}
	}

	//hte��hre
	hte=hts-hst;
	hre=hrs-hsr;
}




/*
�������ܣ�
·����������� ����ģ�Ͳ���
*/
void CP1812_Fun::DeriDiffPar(double hst,double hsr,double htc,double hrc,int n,double *dn,double *hn,
						 double &hstd, double &hsrd)
{
	//����Ժʿ�Ĳ�ͬ: hn->gn
	double d1,d2;
	double d=*(dn+n-1);
	double hobs,aobt,aobr;
	double gt,gr,hstp,hsrp,hi,Hi,Hj,Hk;
	double h1=*hn;
	double hnn=*(hn+n-1);
	
	hobs=-999;
	aobt=-999;
	aobr=-999;
	for (int i=1;i<n-1;i++)
	{
		d1=*(dn+i);
		d2=d-d1;
		hi=*(hn+i);
		
		Hi=hi-(htc*d2+hrc*d1)/d;
		Hj=Hi/d1;
		Hk=Hi/d2;
		
		if (Hi>hobs)//�ж�
		{
			hobs=Hi;
		}	
		if (Hj>aobt)//�ж�
		{
			aobt=Hj;
		}		
		if (Hk>aobr)//�ж�
		{
			aobr=Hk;
		}
		
	}
	
	gt=aobt/(aobt+aobr);
	gr=aobr/(aobt+aobr);
	
	if(hobs>0)
	{
		hstp=hst-hobs*gt;
		hsrp=hsr-hobs*gr;
	}
	else
	{
		hstp=hst;
		hsrp=hsr;
	}
	
	
	hstd=min(h1, hstp);
	hsrd=min(hnn, hsrp);
	
	
	//htesph=htc-hstd;
	//hresph=hrc-hsrd;
	
}




/*
�������ܣ����㲨�����	
*/
double CP1812_Fun::LosDuct(double freq,double p,double b0,double hts,double hrs,int n, double *dn,double *hn, int *climn, 
							  double ae,double w,double tau, int ilt, int ilr, double thetat, double thetar) 
{	
	double dct,dcr;
	GetDctDcr(n,dn,climn, dct,dcr);

    double hst,hsr;
	DeriSmoothEarthsurface(n,dn,hn, hst,hsr);

	double hm,hte,hre;
	DeriDuctPar(n,dn,hn, hts,hrs,ilt,ilr,hst,hsr,hm,hte,hre);

	
	double d=*(dn+n-1);
	double dlt,dlr;
	dlt=*(dn+ilt);
	dlr=d-*(dn+ilr);
	

	//_________________�õ�Af_______________________//
	double Af;
	double thetatpp=thetat-0.1*dlt;
	double thetarpp=thetar-0.1*dlr;

	double Alf=0;//(zr)???�ش���ģ�ԭ���ĳ�������û�������������֪����ô����©���ˣ�
	if (freq<0.5)
	{
		Alf=45.375-137.*freq+92.5*freq*freq; 
	}

	double Ast=0;
	double Asr=0;;
	if (thetatpp>0)
	{
		Ast=20*log10( 1+0.361*thetatpp*sqrt(freq*dlt) )+ 0.264*thetatpp*pow(freq,1./3.);
	} 

	if (thetarpp>0)
	{
		Asr=20*log10( 1+0.361*thetarpp*sqrt(freq*dlr) )+ 0.264*thetarpp*pow(freq,1./3.);
	} 


	double Act=0;
	double Acr=0;
	if ( (w>=0.75) && (dct<=dlt) && (dct<=5.) )
	{
		Act=-3.*exp(-0.25*dct*dct)*(1+tanh(0.07*(50-hts)));
	} 

	if ( (w>=0.75) && (dcr<=dlr) && (dcr<=5.) )
	{
		Acr=-3.*exp(-0.25*dcr*dcr)*(1+tanh(0.07*(50-hrs)));
	} 


	Af=102.45+20*log10(freq)+20*log10(dlt+dlr)+Alf+Ast+Asr+Act+Acr;

	//_________________�õ�Ad(p)_______________________//
	double Adp;
	double gammad=5.e-5*ae*pow(freq,1./3.);

	double thetap,thetatp,thetarp;
	thetatp=min( thetat,0.1*dlt );
	thetarp=min( thetar,0.1*dlr );
	thetap=1.e3*d/ae+thetatp+thetarp;

	double u3;
	double dI=min(d-dlt-dlr, 40.);
	if (hm<=10) 
	{
		u3=1;
	} 
	else
	{
		u3=exp( -4.6e-5*(hm-10.)*(43.+6.*dI) );
	}

	double epsilon=3.5;
	double alpha=-0.6-epsilon*1.e-9*pow(d,3.1)*tau;
	alpha=max(alpha,-3.4);
	double umid=(500./ae) *pow(d,2)/pow(sqrt(hte)+sqrt(hre),2);
	double u2=pow(umid,alpha);
	u2=min(u2,1.);

	double beta=b0*u2*u3;

	//(zr)???���ݵ���������ȷ��ָ��,�������޸ģ���Ϊ��ʱ�����ֵ���Ϊ�������
	double bnum=2.0058-log10(beta);
	double inum;
	if (bnum>0)
	{
		inum=1.012;
	} 
	else
	{
		inum=1;
	}
	double m1=1.076/pow(bnum,inum);

	double m2=-1.e-6*( 9.51-4.8*log10(beta)+0.198*log10(beta)*log10(beta) )*pow(d,1.13);
	double Gamma=m1*exp(m2);

	double Ap=-12.+ (1.2+3.7e-3*d)*log10(p/beta)+12.*pow(p/beta,Gamma);

	Adp=gammad*thetap+Ap;

	//_________________�õ�Lba_______________________//
	double lba;
	lba=Af+Adp;


	
	return lba;
}


/*
�������ܣ���ȡdct/dcr
�����������õõ�������dct/dcr���շ�վ�������ߵľ��룬���Ǻ��ϵľ���, �õ����������Ϣclimn����fortran��һ��
*/
void CP1812_Fun::GetDctDcr(int n, double *dn, int *climn, double &dct, double &dcr)
{	
	for (int i=0;i<n;i++)
	{
		if ( *(climn+i)==1 || *(dn+i)>5 || (i==n-1) ) //�������ȽϽ������
		{
			dct=*(dn+i);
			break;
		} 
	}
	
	for (i=n-1;i>=0;i--)
	{
		if (*(climn+i)==1 || (*(dn+n-1)-*(dn+i))>5 || (i==0) )
		{
			dcr=*(dn+n-1)-*(dn+i);
			break;
		}
	}

	
}


double CP1812_Fun::Stim(double htc,int n, double *dn, double *gn, double ae)//��ԭ�����е�Im50�����޸ĵõ�
{
	double d=*(dn+n-1);
	double d1,d2,gi;
	double nu,numax;
	
	//�������濼��0-n
	numax=-999; //��ʼֵ������С����
	for (int i=1;i<n-1;i++)
	{
		d1=*(dn+i);
		d2=d-d1;
		gi=*(gn+i);
		nu=(gi+500*d1*d2/ae-htc)/d1;
		
		if (nu>numax)//�ж�
		{
			numax=nu;
		}		
	}
	
	return numax;
}

double CP1812_Fun::Vmax(double freq,double h1,double h2, int n, double *dn, double *gn, double ae)
{
	
	double d1,d2,gi;
	double nu,numax;
	double lamda=0.3/freq;
	double d=*(dn+n-1);
	
	//�������濼��0-n
	numax=-999; //��ʼֵ������С����
	for (int i=1;i<n-1;i++)
	{
		d1=*(dn+i);
		d2=d-d1;
		gi=*(gn+i);
		nu=(gi+500*d1*d2/ae-(h1*d2+h2*d1)/d)*sqrt(0.002*d/lamda/d1/d2);
		
		if (nu>numax)//�ж�
		{
			numax=nu;
		}		
	}
	
	return numax;
}


//compute one knife-edge loss
double CP1812_Fun::LosKnifeEdge(double v)
{
	
	double st,jv;
	if (v>-0.78)
	{
		st=sqrt( pow(v-0.1,2)+1 ); 
		jv=6.9+20*log10(st+v-0.1);
	} 
	else
	{
		jv=0;
	}
	
	return jv;
	
}

double CP1812_Fun::Srim(double hrc, int n, double *dn, double *gn, double ae)//��ԭ�����е�Im50�����޸ĵõ�
{
	double d=*(dn+n-1);
	double d1,d2,gi;
	double nu,numax;
	
	//�������濼��0-n
	numax=-999; //��ʼֵ������С����
	for (int i=1;i<n-1;i++)
	{
		d1=*(dn+i);
		d2=d-d1;
		gi=*(gn+i);
		nu=(gi+500*d1*d2/ae-hrc)/d2;
		
		if (nu>numax)//�ж�
		{
			numax=nu;
		}		
	}
	
	return numax;
}

/*
�������ܣ�
4.3.3  First-term part of spherical-Earth diffraction loss��ˮ���½�صĸ������
*/
double CP1812_Fun::LdFirstTermComp(double freq, double d,double htesph,double hresph, double adft,double w,int pol)
{
	double eps1=22.0;
	double sgm1=0.003;
	double eps2=80.0;
	double sgm2=5.0;
	double ldftland=LdFirstTerm(freq, d, htesph, hresph, eps1, sgm1, adft, pol);
	double ldftsea =LdFirstTerm(freq, d, htesph, hresph, eps2, sgm2, adft, pol);
	double ldft=w*ldftsea+(1-w)*ldftland;
	
	return ldft;
	
}


/*
����˵����
4.3.3  First-term part of spherical-Earth diffraction loss; �������;

����˵����
  
���ߺ����ڣ�
	
����˵����
*/
double CP1812_Fun::LdFirstTerm(double freq, double d,double htesph,double hresph,double eps,double sgm, double adft,int pol)
{
	double ldft;
	double k;
	k=0.036*pow(adft*freq,-1./3.)*pow( (pow(eps-1,2)+pow(18*sgm/freq,2)), -1./4.);
	if(pol==1) // pol==1Ϊ��ֱ����
	{
		k=k*pow( eps*eps+pow(18*sgm/freq,2), 1./2.);
	}

	double beta1=1+1.6*pow(k,2)+0.67*pow(k,4);
	double beta2=1+4.5*pow(k,2)+1.53*pow(k,4);
	double beta=beta1/beta2;

	double x=21.88*beta*pow(freq,1./3.)*pow(adft,-2./3.)*d;
	double y=0.9575*beta*pow(freq,2./3.)*pow(adft,-1./3.);
	double y1=y*htesph;
	double y2=y*hresph;

	double fx;
	double b1=beta*y1;
	double b2=beta*y2;
	if(x<1.6)
	{fx=-20*log10(x)-5.6488*pow(x,1.425);}
	else
	{fx=11+10*log10(x)-17.6*x; }

	double gy1,gy2;
	if(b1>2)
	{gy1=17.6*pow(b1-1.1,0.5)-5*log10(b1-1.1)-8;}
	else
	{gy1=20*log10(b1+0.1*pow(b1,3));}

	if(b2>2)
	{gy2=17.6*pow(b2-1.1,0.5)-5*log10(b2-1.1)-8;}
	else
	{gy2=20*log10(b2+0.1*pow(b2,3));}

	ldft=-fx-gy1-gy2;

	return ldft;

}


double CP1812_Fun::LosTropoScatter(double freq, double p, double d, double theta, double N0)
{
	//���������ɢ����ģ��õ���itu�еĹ�ʽ
	
	double lf,lbs;
	lf=25*log10(freq)-2.5*log10(freq/2)*log10(freq/2);
	lbs=190.1+lf+20*log10(d)+0.573*theta-0.15*N0-10.125*pow(log10(50/p), 0.7);
	
	return lbs;
}

double CP1812_Fun::GetFjk(double Theta, double xi, double theta)
{
	//��ֵ����Fj��Fk�ĺ���
	double F=1-0.5*( 1+tanh( 3*xi*(theta-Theta)/Theta ) );
	return F;
	
}


//�ۻ���̬�ֲ����溯��
double CP1812_Fun::I(double x)
{
	
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


double CP1812_Fun::CompLossLoc(double freq, double pL, double h, double R, double lb0p,double lbc)
{
	double sigmaL = (0.024*freq+0.52) * pow(100.,0.28);
	
	double uh;
	if (h<R) //h>=0
	{
		uh=1;
	}
	else if (h<R+10)
	{
		uh=1-(h-R)/10.;
	}
	else
	{
		uh=0;
	}
	
	double sigmaLoc;
	sigmaLoc = uh*sigmaL;
	double Lloc=0;
	double Lb = max(lb0p,lbc+Lloc-I(pL/100.)*sigmaLoc);

	return Lb;
}


/*
����˵�����������ת��Ϊ��ǿ, Ĭ����Ч���书��Ϊ1kW
*/
double CP1812_Fun::L2E(double freq, double lb)
{
	
	double Ep;
	Ep=199.36+20*log10(freq)-lb;
	
	return Ep;
}

//�ۻ���̬�ֲ����溯��__T(x)
double CP1812_Fun::T(double x)
{
	double y=sqrt(-2*log(x));
	return y;
}

double CP1812_Fun::Xi(double x)
{
	//�ۻ���̬�ֲ����溯��__Xi(x)
	
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



/*
����˵����
P1812-6������������ITU��վ�ϵ�һ�£����أ��������࣬���·�������Ϣ��

���������
freq����Ƶ�ʣ�GHz
p����ʱ��ٷֱȣ�p%
pL�����ص�ٷֱȣ�pL%
latMid����·���е��γ�ȣ���
htg��������վ�����߸߶ȣ�agl. m
hrg��������վ�����߸߶ȣ�agl. m
pol����������0:ˮƽ,1:��ֱ
erp������Ч���书�ʣ���kW����λ�ȽϺ��ʣ�(kW)

n����n>=3
dn����dn�ļ����һ���ǵȼ��ģ�ע�ⵥλ��km
hn����m
Rn�������ӵ���߶ȣ�m
climn�����������Ϣ�������P452��һ�£�Ҳ�������ԭ�������еı�����Coastal land A1(1), Inland A2 (2) or Sea B (3);
dN����������1km���������ݶȣ�dN>0����Ҫȡ��ֵ��
N0����


���������
ae������Ч����뾶 (km)
bpath�����Ӿ�·��Ϊtrue,���Ӿ�·��Ϊfalse
thetat�����������ƽ������ (mrad)
theta����·���Ǿ��� (mrad)
thetar�������ն���ƽ������ (mrad)
hm�������δֲڶ� (m)
Lb����������� (dB)	
Ep���� ���ճ�ǿ (dB(uV/m))
errorMes����������Ϣ
*/
bool CP1812_Fun::P1812_6_ITUwb(double freq, double p, double pL, double latMid, double htg, double hrg, int pol, double erp,//�������е�����
							   int n, double *dn, double *hn,double *Rn, int *climn, double dN, double N0, //�߳����ݺͻ�������
							   double &ae, bool &bpath, double &thetat, double &thetar, double &theta, double &hm, double &Lb,double &Ep, CString &errorMes) //�������					  
							 
{
	//������������������(zr)?����������ʱ���ټӣ���������//	
	if (freq<0.03 || freq>6)
	{
		errorMes="Ƶ�ʷ�ΧΪ[0.03, 6]GHz��";
		return false;
	}
	if (p<1 || p>50)
	{
		errorMes="ʱ����ʷ�ΧΪ[1, 50]%��";
		return false;
	}

	if (pL<1 || pL>99)
	{
		errorMes="�ص���ʷ�ΧΪ[1, 99]%��";
		return false;
	}

	if (latMid<-80 || latMid>80)
	{
		errorMes="����·���е�γ�ȷ�ΧΪ[-80, 80]�ȣ�";
		return false;
	}

	if (htg<=0 || htg>3000)
	{
		errorMes="�������ߵĸ߶ȷ�ΧΪ(0, 3000]m��";
		return false;
	}

	if (hrg<=0 || hrg>3000)
	{
		errorMes="�������ߵĸ߶ȷ�ΧΪ(0, 3000]m��";
		return false;
	}

	if (pol!=0 && pol!=1)
	{
		errorMes="���߼�����ȡֵΪ0��1��";
		return false;
	}
	if (n<3)
	{
		errorMes="·�������������Ϊ3����";
		return false;
	}
	

	//������������Ԥ����������//
	double *gn=new double[n];
	*gn = *hn;
	for (int i=1;i<n-1;i++) //����˺ͽ��ն˲��Ӹ��ӵ���߶�
	{
		*(gn+i)=*(hn+i) + *(Rn+i);
	}
	*(gn+n-1) =*(hn+n-1);		 

	
	
	//���������������㡪������//
	dN = fabs(dN);
	double b0,tau,w; 
	ParaInit(latMid,n,dn,climn,dN,ae,b0,tau,w);


	//(zr)������������������㡪������//
	double hts=*hn+htg;
	double hrs=*(hn+n-1)+hrg;
	double htc=hts;	    //(zr)�Ķ���//hts/hrs=htc/hrc 
	double hrc=hrs;
	double d=*(dn+n-1);
	
	//·���������
	int ilt,ilr;
	PathAnalyze(freq,hts,hrs,ae,n,dn,hn, bpath,ilt,ilr,thetat,thetar,theta);
	
    double hst,hsr;
	DeriSmoothEarthsurface(n,dn,hn, hst,hsr);
	
	double hte,hre;
	DeriDuctPar(n,dn,hn, hts,hrs,ilt,ilr,hst,hsr,hm,hte,hre);



	//��������ʱ��仯��������//
	double Lb0p,Lbc;
	CompLossTime(freq,p,pol,htg,hrg,n,dn,hn,gn,climn,N0,ae,b0,tau,w,Lb0p,Lbc);

	//���������ص�仯��������//
	Lb = CompLossLoc(freq,pL,hrg,*(Rn+n-1),Lb0p,Lbc);
		
		
	//����������ġ�>��ǿ��������//
	Ep=L2E(freq,Lb);	
	Ep=Ep+10*log10(erp);


	//ɾ������
	delete []gn;

	return true;
}

/*
����˵�������ݸ߳���������������ȡ���ߵ��������Ϣ

����˵����
climn����radio-climatic zones, 1/A1(coastal land) / 2/A2(inland) / 3/B(sea)���������dct/dcr��dtm/dlm,
ngc����0/2, /0 unknown,1 open/rural,2 water,3 trees/forest,4 suburban, 5 urban,6 dense urban/
  
���ߺ����ڣ�
	
����˵����
2022-5-20���������ƴ���Ŀ��climnΪ���Ϊ����
*/
void CP1812_Fun::RadioClimZone1(int n,double *dn,double *hn,int *ngc,int *climn)
{
	/*
	ȷ��radio-climatic zones,�������dct/dcr��dtm/dlm,1/A1(coastal land) / 2/A2(inland) / 3/B(sea)
	������Ϊ����Ϊ���Ӳ�����ngc/0 unknown,1 open/rural,2 water,3 trees/forest,4 suburban, 5 urban,6 dense urban/
	�����Ӳ����� 2 Ϊˮ�棬�����Ǻ��θ߶�
	
	rczn:�������������

    ˵����(zr)test <=100�ĳ�<100,��������Щ���
	*/

	for (int i=0;i<n;i++)
	{
		*(climn+i)=2; //һ��ʼ��д�� Inland(A2)
	}
	
	
	//ȷ�����������,��չ��
	
	for (i=0;i<n;i++) 
	{
		
		if (*(ngc+i)==2) //�ж��Ƿ�Ϊ����
		{
			*(climn+i)=3;
			
			if ( *(ngc+i-1)!=2 ) //��ǰ��
			{
				int j=1;
				while ( (i-j)>=0 && *(climn+i-j)==2 && *(hn+i-j)<100 && (*(dn+i)-*(dn+i-j))<50  )//(zr)test
				{
					*(climn+i-j)=1;
					j++;
				}				
			} 
			
			else if (*(ngc+i+1)!=2) //�����
			{
				int j=1;
				while ( (i+j)<n && *(climn+i+j)==2 && *(hn+i+j)<100 && *(ngc+i+j)!=2 && (*(dn+i+j)-*(dn+i))<50 )//(zr)test
				{
					*(climn+i+j)=1;
					j++;
				}
			}
		} 	
		
	}
	
}