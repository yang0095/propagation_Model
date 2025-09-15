#include"hazeatt.h"
#include<cmath>
/**
 * 计算霾衰减
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
 * @param haze_type 气象条件（霾）
 * @return 霾衰减
 */
double HazeAtt::haze_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                  Receive_polar receive_polar,double dia,double antenna_effeciency, Haze_type haze_type){
    //检测输入参数
    if(parameter_rational(f,PP,satlon,satlat,satalt,send_polar, lon, lat, alt,receive_polar, dia, antenna_effeciency)==false){
        exit(0);
    }

    double haze;
    if(0==haze_type){
        haze=0.1;
    }else if(1==haze_type){
        haze=0.2;
    }else if(2==haze_type){
        haze=0.4;
    }else{
        haze=0.6;
    }
    haze*=2;
    return std::round(haze * 10) / 10;
}

//std::vector<double> HazeAtt::FUN_m_haze_att(double f, double satlon, double satlat, double satalt, Send_polar send_polar, double ag, double apd, double antenna_effeciency, double lon, double lat, double alt, double Rcs, Haze_type haze_type, std::vector<std::vector<double> > &LatlonArea)
//{
//    std::vector<double> v_result;
//    if(LatlonArea.empty())
//        return v_result;

//    for(std::vector<std::vector<double>>::iterator it=LatlonArea.begin();it!=LatlonArea.end();it++ )
//    {
//        double d_HazeAtt=FUN_haze_att( f,  satlon,  satlat,  satalt,  send_polar,  ag,  apd,  antenna_effeciency,  lon,  lat,  alt,  Rcs,  haze_type,  (*it)[0], (*it)[1], (*it)[3], (*it)[4], (*it)[5]);
//        v_result.push_back(d_HazeAtt);
//    }
//    return v_result;
//}

double HazeAtt::FUN_haze_att(double f,double satlat,double satlon, double satalt,
                             Send_polar send_polar, double antenna_effeciency,
                             double lat, double lon, double alt,int haze_area_count,
                             vector<vector<double>>haze_area, vector<Haze_type>hR)
{
    RainAtt myrainatt;
    CcalculateElevation mycalele;
    double ele= mycalele.calculateElevation(satlat, satlon,satalt,lat,lon,alt);
    double haze_att_count=0.0;

    for(int i=0;i<haze_area_count;i++){
        if(hR[i]==Haze_type::Slight){
            continue;
        }
        //检测输入参数
        if(parameter_rational1(f,satlon,satlat,satalt,send_polar, lon, lat, alt, antenna_effeciency)==false){
            return 0;
        }
        try {
                // 计算有效路径长度

                double r = myrainatt.calculateRainPathLength(
                    satlat, satlon, satalt, lat, lon, alt,
                    haze_area[i][0], haze_area[i][1], haze_area[i][3], haze_area[i][4],
                    haze_area[i][5], ele );

                std::cout << "霾区内有效路径长度: " << r << " 米\n";

                r/=1000;//(km)
                // 计算衰减
                double haze=specificatten(f,hR[i]);
                std::cout << "霾衰减率: " << haze ;
                haze_att_count+= haze * r;

            } catch (const std::exception& e) {
                std::cerr << "错误: " << e.what() << std::endl;
                return 1;
            }

    }
    cout<<"霾总衰减"<<haze_att_count<<endl;
    return haze_att_count;



}



bool HazeAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                Receive_polar receive_polar,double dia,double antenna_effeciency){

    bool flag=true;
    f=f/pow(10,9);
    if(f<12||f>18){
        std::cout<<"haze_att:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"haze_att:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"haze_att:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"haze_att:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"haze_att:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"haze_att:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(dia<0||dia>5){
        std::cout<<"haze_att:天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
        exit(0);
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"haze_att:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }

    return flag;
}

bool HazeAtt::parameter_rational1(double f, double satlon, double satlat, double satalt, Send_polar send_polar, double lon, double lat, double alt, double antenna_effeciency)
{
    bool flag=true;
//    f=f/pow(10,9);
    if(f<12||f>38){
        std::cout<<"haze_att:频率="<<f<<"(GHZ) 超出限定范围（12~18，32~38）!"<<endl;
        flag=false;
    }else if(f>18 && f<32)
    {
        std::cout<<"haze_att:频率="<<f<<"(GHZ) 超出限定范围（12~18，32~38）!"<<endl;
        flag=false;
    }
    if(satlon<-180||satlon>180){
        std::cout<<"haze_att:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"haze_att:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"haze_att:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"haze_att:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"haze_att:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"haze_att:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }

    return flag;
}

double HazeAtt::specificatten(double f, Haze_type haze_type)
{
    double haze=0.0;  //（zr）插值计算
    //zy回复：这里使用线性插值，根据频点表进行重新计算,此处存疑，无返回值？
    //先计算Ku频段范围12-18Ghz,然后是Ka频段32-38Ghz
    //（zr）-1表示没有
    if(f>=12.0 && f<13){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000003;
        }else if(1==haze_type){
            haze=(f-12)*(0.000235-0.000217)+0.000217;
        }else if(2==haze_type){
            haze=(f-12)*(0.004888-0.004512)+0.004512;
        }
    }else if(f>=13 &&f<14){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000003;
        }else if(1==haze_type){
            haze=(f-13)*(0.000253-0.000235)+0.000235;
        }else if(2==haze_type){
            haze=(f-13)*(0.005264-0.004888)+0.004888;
        }
    }else if(f>=14 &&f<15){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000003;
        }else if(1==haze_type){
            haze=(f-14)*(0.000271-0.000253)+0.000253;
        }else if(2==haze_type){
            haze=(f-14)*(0.005640-0.005264)+0.005264;
        }
    }else if(f>=15 &&f<16){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000003;
        }else if(1==haze_type){
            haze=(f-15)*(0.000289-0.000271)+0.000271;
        }else if(2==haze_type){
            haze=(f-15)*(0.006016-0.005640)+0.005640;
        }
    }else if(f>=16 &&f<17){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=(f-16)*(0.000004-0.000003)+0.000003;
        }else if(1==haze_type){
            haze=(f-16)*(0.000307-0.000289)+0.000289;
        }else if(2==haze_type){
            haze=(f-16)*(0.006393-0.006016)+0.006016;
        }
    }else if(f>=17 &&f<18){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000004;
        }else if(1==haze_type){
            haze=(f-17)*(0.000325-0.000307)+0.000307;
        }else if(2==haze_type){
            haze=(f-17)*(0.006769-0.006393)+0.006393;
        }
    }else if(f>=32 &&f<33){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000007;
        }else if(1==haze_type){
            haze=(f-32)*(0.000596-0.000578)+0.000578;
        }else if(2==haze_type){
            haze=(f-32)*(0.012033-0.012409)+0.012033;
        }
    }else if(f>=33 &&f<34){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000007;
        }else if(1==haze_type){
            haze=(f-33)*(0.000615-0.000596)+0.000596;
        }else if(2==haze_type){
            haze=(f-33)*(0.012785-0.012409)+0.012409;
        }
    }else if(f>=34 &&f<35){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=(f-34)*(0.000008-0.000007)+0.000007;
        }else if(1==haze_type){
            haze=(f-34)*(0.000633-0.000615)+0.000615;
        }else if(2==haze_type){
            haze=(f-34)*(0.013161-0.012785)+0.012785;
        }
    }else if(f>=35 &&f<36){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000008;
        }else if(1==haze_type){
            haze=(f-35)*(0.000651-0.000633)+0.000633;
        }else if(2==haze_type){
            haze=(f-35)*(0.013538-0.013161)+0.013161;
        }
    }else if(f>=36 &&f<37){
        if(-1==haze_type){
            haze=0.0;
        }else if(0==haze_type){
            haze=0.000008;
        }else if(1==haze_type){
            haze=(f-36)*(0.000669-0.000651)+0.000651;
        }else if(2==haze_type){
            haze=(f-36)*(0.013914-0.013538)+0.013538;
        }
    }
    else if(f>=37 &&f<=38){
            if(-1==haze_type){
                haze=0.0;
            }else if(0==haze_type){
                haze=0.000008;
            }else if(1==haze_type){
                haze=(f-37)*(0.000687-0.000669)+0.000669;
            }else if(2==haze_type){
                haze=(f-37)*(0.014290-0.013914)+0.013914;
            }
        }
    return haze;


//    if(f==12.0)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000003;
//        }else if(2==haze_type){
//            haze=0.000217;
//        }else{
//            haze=0.004512;
//        }
//    }else if(f==13)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000003;
//        }else if(2==haze_type){
//            haze=0.000235;
//        }else{
//            haze=0.004888;
//        }
//    }else if(f==14)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000003;
//        }else if(2==haze_type){
//            haze=0.000253;
//        }else{
//            haze=0.005264;
//        }
//    }else if(f==15)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000003;
//        }else if(2==haze_type){
//            haze=0.000271;
//        }else{
//            haze=0.005640;
//        }
//    }else if(f==16)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000003;
//        }else if(2==haze_type){
//            haze=0.000289;
//        }else{
//            haze=0.006061;
//        }
//    }else if(f==17)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000004;
//        }else if(2==haze_type){
//            haze=0.000307;
//        }else{
//            haze=0.006393;
//        }
//    }else if(f==18)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000004;
//        }else if(2==haze_type){
//            haze=0.000325;
//        }else{
//            haze=0.006769;
//        }
//    }else if(f==32)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000007;
//        }else if(2==haze_type){
//            haze=0.000578;
//        }else{
//            haze=0.012033;
//        }
//    }else if(f==33)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000007;
//        }else if(2==haze_type){
//            haze=0.000596;
//        }else{
//            haze=0.012409;
//        }
//    }else if(f==34)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000007;
//        }else if(2==haze_type){
//            haze=0.000615;
//        }else{
//            haze=0.012785;
//        }
//    }else if(f==35)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000008;
//        }else if(2==haze_type){
//            haze=0.000633;
//        }else{
//            haze=0.013161;
//        }
//    }else if(f==36)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000008;
//        }else if(2==haze_type){
//            haze=0.000651;
//        }else{
//            haze=0.013538;
//        }
//    }else if(f==37)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000008;
//        }else if(2==haze_type){
//            haze=0.000669;
//        }else{
//            haze=0.013914;
//        }
//    }else if(f==38)
//    {
//        if(0==haze_type){
//            haze=0.0;
//        }else if(1==haze_type){
//            haze=0.000008;
//        }else if(2==haze_type){
//            haze=0.000687;
//        }else{
//            haze=0.014290;
//        }
//    }

}


