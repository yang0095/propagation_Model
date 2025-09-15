#include"fogatt.h"

/**
 * 计算雾衰减
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
 *
 * @param fog_meteor_condtion //气象条件（雾）等级
 * @param start_distance   气象条件（雾）起始距离
 * @param end_distance 气象条件（雾）终止距离
 * @param L 液态水含量(mm) 精度0.10.1
 * @param water_tmp 雾中液态水温度（K）精度0.1
 * @return 雾衰减
 */
double FogAtt::fog_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
               Receive_polar receive_polar,double dia,double antenna_effeciency, Fog_meteor_condtion meteor_condtion,double start_distance,double end_distance,double L,double water_tmp){
    //判断参数合理性
    if(parameter_rational(f, PP, satlon,  satlat, satalt,  send_polar, lon,  lat, alt, receive_polar, dia, antenna_effeciency,  meteor_condtion, start_distance, end_distance,L, water_tmp)==false){
        exit(0);
    }

  if(0==meteor_condtion){
        L=std::max(L,0.05);
    }else{
        L=std::max(L,0.5);
    }
    CloudAtt cloud;
    //得到云衰减
    Cloud_type cloud_type=static_cast<Cloud_type>(0);
    double cloud_att=cloud.cloud_att(f,PP,satlon,satlat,satalt,send_polar,lon,lat,alt,receive_polar,dia,antenna_effeciency,cloud_type,L,water_tmp);
    double result=cloud_att*L*abs(start_distance-end_distance)/1000;
    result*=2;
    return std::round(result * 10) / 10;
}

double FogAtt::FUN_fog_att(double f,double temp,
                           double lat1, double lon1, double h1,
                           double lat2, double lon2, double h2,
                           int fog_area_count,vector<vector<double>>fog_area, vector<double > fR)
{
    //（zr）water_tmp=地面温度
    double water_tmp=temp;//zy没有水温输入参数，直接等于地面温度

    RainAtt myrain;
    CloudAtt mycloud;
    CcalculateElevation mycalelev;
    double elev_d;

    // 计算仰角
    elev_d =mycalelev.calculateElevation(lat1,lon1,h1,lat2,lon2,h2);
    cout<<"计算角度为"<<elev_d<<endl;
    double att=0.0;
    for(int i=0;i<fog_area_count;i++){
        if(parameter_rational1(f,fR[i],water_tmp,lat1, lon1,lat2, lon2,
                               fog_area[i][0], fog_area[i][1], fog_area[i][3], fog_area[i][4])==false){
            exit(0);
        }

        try {
                // 计算有效路径长度
                double r_count = myrain.calculateRainPathLength(
                    lat1, lon1, h1, lat2, lon2, h2,
                    fog_area[i][0], fog_area[i][1], fog_area[i][3], fog_area[i][4],
                    fog_area[i][5], elev_d);

                std::cout << "雾区内有效路径长度: " << r_count << " 米\n";
                double atten_km=mycloud.specificatten(f,  fR[i], water_tmp);
                std::cout << "雾衰减率: " << atten_km ;
                // 计算衰减
                att+=r_count*atten_km/1000;




            } catch (const std::exception& e) {
                std::cerr << "错误: " << e.what() << std::endl;
                return 1;
            }



    }
    std::cout << "雾衰减: " << att<<endl ;
    return att;



        return 0;
}



bool FogAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                Receive_polar receive_polar,double dia,double antenna_effeciency, Fog_meteor_condtion meteor_condtion,double start_distance,double end_distance,double L,double water_tmp){

    bool flag=true;
    f=f/pow(10,9);
    if(f<12||f>18){
        std::cout<<"fogatt:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"fogatt:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"fogatt:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"fogatt:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"fogatt:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"fogatt:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(dia<0||dia>5){
        std::cout<<"fogatt:天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
        exit(0);
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"fogatt:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }
    if(L<0||L>2){
        std::cout<<"fogatt:液态水含量="<<L<<"超出限定范围（0~2）!"<<endl;
        flag=false;
    }

    if(water_tmp<193.5||water_tmp>323.5){
        std::cout<<"fogatt:云中液态水温度="<<water_tmp<<"超出限定范围（193.5~323.5）!"<<endl;
        flag=false;
    }

    return flag;
}

bool FogAtt::parameter_rational1(double f, double L, double water_tmp, double lat1, double lon1, double lat2, double lon2, double lata, double lona, double latb, double lonb)
{
    bool flag=true;
    if(f<12||f>18){
        std::cout<<"FogAtt:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }
    if(L<0||L>10){
        std::cout<<"FogAtt:液态水含量="<<L<<"超出限定范围（0~10）!"<<endl;
        flag=false;
    }

    if(water_tmp<193.5||water_tmp>323.5){
        std::cout<<"FogAtt:雾中液态水温度="<<water_tmp<<"超出限定范围（193.5~323.5）!"<<endl;
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
