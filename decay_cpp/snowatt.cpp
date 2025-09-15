#include"snowatt.h"


/**
 * 计算雪衰减
 * @param f 频率f（GHz） 精度1
 * @param PP 时间概率（%） 精度0.001
 * @param satlon 卫星经度（度）精度10-6
 * @param satlat 卫星纬度（度）精度10-6
 * @param satalt 飞行器高度（m）精度0.1
 * @param send_polar  信源发射天线极化方式
 * @param lon 站点经度（度）精度10-6
 * @param lat 站点纬度（度）精度10-6
 * @param alt 站点高度（单位：m) 精度0.1
 * @param receive_polar  信宿接收天线极化方式
 * @param dia 直径（m） 精度0.1
 * @param antenna_effeciency 天线效率 精度0.01
 * @param snow_type 气象条件（雪）降雪量
 * @param start_distance   气象条件（雪）起始距离
 * @param end_distance 气象条件（雪）终止距离
 * @return 雪衰减（dB）
 */
double SnowAtt::snow_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                  Receive_polar receive_polar,double dia,double antenna_effeciency, Snow_type snow_type,double start_distance,double end_distance){
    //判断参数合理性
    if(parameter_rational(f, PP, satlon, satlat, satalt, send_polar,lon, lat,alt,receive_polar, dia, antenna_effeciency)==false){
        exit(0);
    }

    f=f/pow(10,9);
    if(f>20.0){
        std::cout<<"雪衰减:输入频率高于20GHz，数据有误！"<<std::endl;
        exit(1);
    }

    double lamda = 30/f; //%波长，cm
    double r;
    if(0==static_cast<int>(snow_type)){
        r=0.1;
    }else if(1==static_cast<int>(snow_type)){
        r=0.5;
    }else{
        r=1;
    }

    double alpha = 0.00349*pow(r,1.6)/pow(lamda,4) + 0.00224*r/lamda;
    double distance=abs(end_distance-start_distance)/1000;//m->km
    double result=alpha*distance;
    result*=2;
    return std::round(result * 10) / 10;
}

double SnowAtt::specificatten(double f, Snow_type snow_type)
{
//    f=f/pow(10,9);
    if(f>20.0){
        std::cout<<"雪衰减:输入频率高于20GHz，数据有误！"<<std::endl;
        exit(1);
    }

    double lamda = 30/f; //%波长，cm
    double r;//降水率
    if(0==static_cast<int>(snow_type)){
        r=0.1;
    }else if(1==static_cast<int>(snow_type)){
        r=0.5;
    }else{
        r=1;
    }

    double alpha = 0.00349*pow(r,1.6)/pow(lamda,4) + 0.00224*r/lamda;
    return alpha;
}

double SnowAtt::FUN_snow_att(double f, double lat1, double lon1, double h1,
                             double lat2, double lon2, double h2,
                             int snow_area_count,vector<vector<double>>snow_area,vector<Snow_type>sR)
{
    RainAtt myrainatt;
    CcalculateElevation mycalele;
    double elev_d= mycalele.calculateElevation(lat1, lon1,h1,lat2,lon2,h2);
    double snow_att_count=0.0;
    for(int i=0;i<snow_area_count;i++){
        if(parameter_rational1(f,lat1, lon1,lat2, lon2, snow_area[i][0], snow_area[i][1], snow_area[i][3], snow_area[i][4])==false){
            exit(0);
        }
        if(sR[i]==Snow_type::None){continue;}
        try {
                // 计算有效路径长度
                double r = myrainatt.calculateRainPathLength(
                    lat1, lon1, h1, lat2, lon2, h2,
                    snow_area[i][0], snow_area[i][1], snow_area[i][3], snow_area[i][4],
                    snow_area[i][5], elev_d);

                std::cout << "雪区内有效路径长度: " << r << " 米\n";

                // 计算衰减
                double atten_km=specificatten(f, sR[i]);
                std::cout << "雪衰减率: " << atten_km ;
                snow_att_count+= atten_km * r/1000;

            } catch (const std::exception& e) {
                std::cerr << "错误: " << e.what() << std::endl;
                return 1;
            }
    }
    cout<<"雪区总衰减"<<snow_att_count<<endl;
        return snow_att_count;

}



bool SnowAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                Receive_polar receive_polar,double dia,double antenna_effeciency){

    bool flag=true;
    f=f/pow(10,9);
    if(f<12||f>18){
        std::cout<<"snow_att:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"snow_att:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"snow_att:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"snow_att:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"snow_att:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"snow_att:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(dia<0||dia>5){
        std::cout<<"snow_att:天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
        exit(0);
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"snow_att:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }

    return flag;
}

bool SnowAtt::parameter_rational1(double f, double lat1, double lon1, double lat2, double lon2, double lata, double lona, double latb, double lonb)
{
    bool flag=true;
    if(f<12||f>18){
        std::cout<<"FogAtt:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }
    if(lon1<-180||lon1>180){
        std::cout<<"FogAtt:位置1经度="<<lon1<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat1<-90||lat1>90){
        std::cout<<"FogAtt:位置1纬度="<<lat1<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(lon2<-180||lon2>180){
        std::cout<<"FogAtt:位置2经度="<<lon2<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat2<-90||lat2>90){
        std::cout<<"FogAtt:位置2纬度="<<lat2<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(lona<-180||lona>180){
        std::cout<<"FogAtt:左下经度="<<lona<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lata<-90||lata>90){
        std::cout<<"FogAtt:左下纬度="<<lata<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(lonb<-180||lonb>180){
        std::cout<<"FogAtt:右上经度="<<lonb<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(latb<-90||latb>90){
        std::cout<<"FogAtt:右上纬度="<<latb<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    return flag;
}

