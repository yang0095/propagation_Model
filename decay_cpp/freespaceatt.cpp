#include"freespaceatt.h"

//计算球面坐标距离，alt1单位为km
double FreeSpaceAtt::haversine(double lat1,double lon1,double alt1,double lat2,double lon2,double alt2){
    double R=6371.0;//地球半径，单位：公里
    lat1=lat1*M_PI/180.0;
    lon1=lon1*M_PI/180.0;
    lat2=lat2*M_PI/180.0;
    lon2=lon2*M_PI/180.0;
    //哈佛幸公式
    double dlat=lat2-lat1;
    double dlon=lon2-lon1;
    double a=sin(dlat/2)*sin(dlat/2)+cos(lat1)*cos(lat2)*sin(dlon/2)*sin(dlon/2);
    double c=2*atan2(sqrt(a),sqrt(1-a));
    double horizontalDistance=R*c;
    double verticalDistance=alt2-alt1;
    return sqrt(horizontalDistance*horizontalDistance+verticalDistance*verticalDistance);
}


/**
 * 计算自由空间传播衰减
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
 * @return 自由空间传播衰减
 */
double FreeSpaceAtt::freespace_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                  Receive_polar receive_polar,double dia,double antenna_effeciency,double radar,double gain){

    //判断参数输入合理性
    if(parameter_rational(f, PP, satlon,  satlat, satalt,  send_polar, lon,  lat, alt, receive_polar, dia, antenna_effeciency)==false){
        exit(0);
    }
    //利用哈佛幸公式计算距离
   double x=haversine(lat,lon,alt/1000,satlat,satlon,satalt/1000);
    //std::cout<<"距离："<<x<<"km"<<endl;
    f=f*pow(10,3);//频率hz->MHz
    cout<<"freeSpace Att Frequence is"<<f<<endl;
   double L=103.4+20*log10(f)+40*log10(x)-2*gain-radar;
   //New  L=103.4+20*log10(f)+40*log10(x)-双倍天线增益-雷达散射截面
   return std::round(L * 10) / 10;
}

double FreeSpaceAtt::FUN_freespace_att(double f, double PP, double satlon, double satlat, double satalt, Send_polar send_polar, double lon, double lat, double alt, Receive_polar receive_polar, double dia, double antenna_effeciency, double sigma)
{
    //判断参数输入合理性
    if(parameter_rational(f, PP, satlon,  satlat, satalt,  send_polar, lon,  lat, alt, receive_polar, dia, antenna_effeciency)==false){
        exit(0);
    }
    //利用哈佛幸公式计算距离
   double x=haversine(lat,lon,alt/1000,satlat,satlon,satalt/1000);
    //std::cout<<"距离："<<x<<"km"<<endl;
    f=f*pow(10,3);//频率hz->MHz
   double L=103.4+20*log10(f)+40*log10(x)-10*log10(sigma);
   return L;
}

bool FreeSpaceAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                      Receive_polar receive_polar,double dia,double antenna_effeciency){

    bool flag=true;

    if(f<12||f>18){
        std::cout<<"FreeSpaceAtt:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"FreeSpaceAtt:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"FreeSpaceAtt:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"FreeSpaceAtt:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"FreeSpaceAtt:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"FreeSpaceAtt:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(dia<0||dia>5){
        std::cout<<"FreeSpaceAtt:天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
        exit(0);
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"FreeSpaceAtt:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }

    return flag;
}

