#include "sciatt.h"


//根据经度，维度查找文件，返回对应n_wet文件中的数据
double SciAtt::findNWet(const std::string& longitudeFile, const std::string& latitudeFile, const std::string& nWetFile, double inputLongitude, double inputLatitude) {
    std::ifstream lonFile(longitudeFile);
    std::ifstream latFile(latitudeFile);
    std::ifstream wetFile(nWetFile);

    if (!lonFile || !latFile || !wetFile) {
//        std::cerr << "无法打开其中一个文件。" << std::endl;
        return -1; // 返回-1表示错误
    }


    std::string lonLine, latLine, nWetLine;
    // 显式声明lambda表达式的类型为std::function
    std::function<bool(const pair<double, double>&, const pair<double, double>&)> cmp =
        [](const pair<double, double>& p1, const pair<double, double>& p2) {
            return p1.first < p2.first; // 大根堆：第一个元素大的优先
        };

    // 创建一个大根堆，存储pair<double, double>，使用显式声明的lambda表达式,设置大小为2
    priority_queue<pair<double, double>, vector<pair<double, double>>, std::function<bool(const pair<double, double>&, const pair<double, double>&)>> maxHeap(cmp);

    //逐行读取数据，将每组数据和输入的经度纬度求曼哈顿距离，保持最小的两个在大根堆中
    while (std::getline(lonFile, lonLine) && std::getline(latFile, latLine) && std::getline(wetFile, nWetLine)) {

        std::istringstream lonStream(lonLine);
        std::istringstream latStream(latLine);
        std::istringstream nWetStream(nWetLine);

        double longitude, latitude, n_wet;

        // 逐个读取行内的多个值
        while (lonStream >> longitude && latStream >> latitude && nWetStream >> n_wet) {
            if(longitude==inputLongitude && latitude==inputLatitude){ //如果精确相等，直接返回对应的n_wet值
                return n_wet;
            }
            double distance=abs(longitude-inputLongitude)+abs(latitude-inputLatitude);
            pair<double,double>pir=make_pair(distance,n_wet);
            if(maxHeap.size()<2){
                maxHeap.push(pir);
            }else if(distance<maxHeap.top().first){
                maxHeap.pop();
                maxHeap.push(pir);
            }
        }
    }
    double arv_n_wet=0.0;
    while(!maxHeap.empty()){
        arv_n_wet+=maxHeap.top().second;
        maxHeap.pop();
    }
    arv_n_wet/=2;
    return arv_n_wet;
}

//格式化文件名
std::string SciAtt::getFileName(double probability) {
    // 将概率转换为 1 到 99 之间的整数
    int probabilityInt = static_cast<int>(std::round(probability * 100));

    // 限制范围在 1 到 99 之间
    if (probabilityInt < 1) {
        probabilityInt = 1;
    } else if (probabilityInt > 99) {
        probabilityInt = 99;
    }

    // 格式化文件名
    std::ostringstream fileNameStream;
    fileNameStream << "NWET_Annual_" << std::setw(2) << std::setfill('0') << probabilityInt << ".TXT";
    return fileNameStream.str();
}



/**
 * 计算闪烁衰减
 * @param freq_GHz 频率f（GHz） 精度1
 * @param PP 时间概率（%） 精度0.001
 * @param lonSatellite 卫星经度（度）精度10-6
 * @param latSatellite 卫星纬度（度）精度10-6
 * @param satalt 飞行器高度（m）精度0.1
 * @param send_polar  信源发射天线极化方式
 * @param lonStation 站点经度（度）精度10-6
 * @param latStation 站点纬度（度）精度10-6
 * @param altStation 站点高度（单位：m) 精度0.1
 * @param receive_polar  信宿接收天线极化方式
 * @param D 直径（m） 精度0.1
 * @param eita 天线效率 精度0.01
 * @return
 */
double SciAtt::sci_att( double f,double PP,double lonSatellite,double latSatellite,double satalt,
                 Send_polar send_polar,double lonStation, double latStation,double altStation, Receive_polar receive_polar,double D,double eita){
    if(parameter_rational( f, PP, lonSatellite, latSatellite, satalt,send_polar, lonStation, latStation, altStation, receive_polar, D, eita)==false){
        exit(0);
    }

      //计算仰角cita
     Position station;
     station.setLat(latStation);
     station.setLon(lonStation);
     station.setAlt(altStation);

     Position satellite;
     satellite.setLat(latSatellite);
     satellite.setLon(lonSatellite);
     satellite.setAlt(satalt);

     Location stationLoc = blhToXyz(station);
     Location satelliteLoc = blhToXyz(satellite);

     double cita = computeAngle(satelliteLoc, stationLoc, station);

     if(cita<5||cita>80){ //判断仰角合理性
         std::cout<<"sciatt:仰角="<<cita<<"不符合范围（5~80）！"<<endl;
         exit(0);
     }

    //获取N_wet文件名
    string nWetFile="NWET_Annual_50.TXT";
    double inputLongitude=(lonSatellite+ lonStation)/2;//求经度平均值
    double inputLatitude=(latSatellite+latStation)/2;//求维度平均值
    string longitudeFile="LON_N.TXT";
    string latitudeFile="LAT_N.TXT";
    ReConfig cfg;
    string src=cfg.getDirPath("sciatt");
    nWetFile=src+nWetFile;
    longitudeFile=src+longitudeFile;
    latitudeFile=src+latitudeFile;

    double  N_wet=findNWet(longitudeFile,latitudeFile,nWetFile,inputLongitude,inputLatitude);

    double a, b, c, d, EF, e_s, e_w, cigema_ref, h_l=1000.0, L, D_eff, A_s, cigema;

    //------------setp-1.计算饱和水汽压；----------------------------------------
    //水；
    // if (type == 1) {
    //     a = 6.1121;, lona, latb, lonb,
    //     b = 18.678;
    //     c = 257.14;
    //     d = 234.5;
    //     EF = 1 + 1e-4 * (7.2 + P * (0.0032 + 5.9e-7 * t * t));
    //     // 冰
    // } else {
    //     a = 6.1115;
    //     b = 23.036;
    //     c = 279.82;
    //     d = 333.7;
    //     EF = 1 + 1e-4 * (2.2 + P * (0.0382 + 6.4e-6 * t * t));
    // }
    // e_s = EF * a * exp(t * (b - t / d) / (t + c));
    //  e_w = rh * e_s / 100; //水汽压强（hPa）；
    //        N_wet = 72 * e_w / T + 3.732e5 * e_w / Math.pow(T, 2); // 折射指数湿项；
    // N_wet = 3.732e5 * e_w / pow(T, 2); // 折射指数湿项；


    // -----------------------3.计算信号幅度的标准偏差；--------------------------
    cigema_ref = 3.6e-3 + 1e-4 * N_wet;
    //---------------------4.计算等效路径长度；----------------------------------
    //（zr）简单处理，等效高度，hl，等效路径长度L=hl/sin(cita) ;PP=50;
    //zy重构部分
    //zy1.闪烁区域范围1km
    //zy2.讨论三种情况：1)二者都大于1Km,则视为等效高度为0
    //zy             2)二者都小于1km，则等效高度L=hl/sin(cita) 其中hl为二者高度差
    //zy             3)一大一小，则计算在1km范围内的等效高度
    double sci_length=1000.0;
    if(satalt>=sci_length &&altStation>=sci_length){
        L=0.0;
    }else if(satalt<sci_length && altStation<sci_length){
        if(satalt>=altStation){
            L=(satalt-altStation)/sin(cita);
        }
    }else if(altStation<sci_length){
        L=(sci_length-altStation)/sin(cita);
    }else{
        L=(sci_length-satalt)/sin(cita);
    }



    //h_l = 2000.0;//修改
    //double h_2 = 1000.0;//修改
    //L = 0.25*abs(h_l-h_2)/(sqrt(sin(cita*M_PI/180)*sin(cita*M_PI/180)+2.35e-4)+sin(cita*M_PI/180));//修改
    //h_l = 1000;    //湍流层高度（m）；
    //L = 2 * h_l / (sqrt(sin(cita * M_PI / 180) * sin(cita * M_PI / 180) + 2.35e-4) + sin(cita * M_PI / 180));

    //---------------------5.估计等效天线尺寸；----------------------------------
    D_eff = sqrt(eita) * D;
    //---------------------6.计算天线平均因子；----------------------------------
    double x = 1.22 * pow(D_eff, 2) * (f / L);
    double A = 3.86 * pow(x * x + 1, 11.0 / 12.0) * sin(11.0 / 6.0 * atan(1.0 / x));
    double B = 7.08 * pow(x, 5.0 / 6.0);
    double g = sqrt(A - B);
    if (A < B || g <= 0) {
        A_s = 0;
    } else {
        // --------------7.计算特定时期和路径的信号标准偏差；---------------------
        cigema = cigema_ref * pow(f, 7.0 / 12.0) * g / pow(sin(cita * M_PI / 180), 1.2);
        // 8.计算时间概率因子；
        a = -0.061 * pow(log10(PP), 3) + 0.072 * pow(log10(PP), 2) - 1.71 * log10(PP) + 3.0;
        // 9.计算p%闪烁衰落深度；
        A_s = a * cigema;
    }
    RainAtt myrainatt;
    double r = myrainatt.calculateRainPathLength(
        latSatellite, lonSatellite, satalt, latStation, lonStation, altStation,
                latSatellite, lonSatellite, latStation, lonStation,
                h_l, cita);
    double Tk =273.15;//开尔文温度
    double waterPressure=Tk/216.7;//水气压
    double wet =72*waterPressure/Tk+3.75*1e5*waterPressure/pow(Tk,2);
    double exponent=7.0/12.0;
    double exponent1=11.0/12.0;
    double exponent2=1.0/3.0;
    wet=0.1*pow(f,exponent)*pow(r/1000,exponent1)/pow(cita,exponent2);
     return wet;
}

double SciAtt::FUN_sci_att(double lat1, double lon1, double h1, double lat2, double lon2, double h2, double lata, double lona, double latb, double lonb, double temperature, double elev_d, double rou, double hR, double PP, double Pressure)
{
   //(zr)没做完
    RainAtt myrainatt;
    if(parameter_rational1(lat1, lon1,lat2, lon2, lata, lona, latb, lonb)==false){
        exit(0);
    }
    try {
            // 计算有效路径长度
            double r = myrainatt.calculateRainPathLength(
                lat1, lon1, h1, lat2, lon2, h2,
                lata, lona, latb, lonb,
                hR, elev_d);

            std::cout << "闪烁区内有效路径长度: " << r << " 米\n";

            // 计算衰减
            double Tk = temperature+273.15;//开尔文温度
            double waterPressure= rou*Tk/216.7;//水气压
            double N_wet =72*waterPressure/Tk+3.75*1e5*waterPressure/pow(Tk,2);

            std::cout << "闪烁衰减: " << N_wet ;
            return N_wet;

        } catch (const std::exception& e) {
            std::cerr << "错误: " << e.what() << std::endl;
            return 1;
        }
        return 0;

}



bool SciAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                 Receive_polar receive_polar,double dia,double antenna_effeciency){

    bool flag=true;

    if(f<12||f>18){
        std::cout<<"sci_att:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"sci_att:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"sci_att:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"sci_att:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"sci_att:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"sci_att:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(dia<0||dia>5){
        std::cout<<"sci_att:天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
        exit(0);
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"sci_att:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }

    return flag;
}

bool SciAtt::parameter_rational1(double lat1, double lon1, double lat2, double lon2, double lata, double lona, double latb, double lonb)
{
    bool flag=true;
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



Location SciAtt::blhToXyz(Position p) {
    double earthRadius = 6378137.0;
    double f = 1.0 / 298.257223563;
    double lat = p.lat * M_PI / 180;
    double lon = p.lon * M_PI / 180;
    double alt = p.alt;

    double e2 = 2 * f - f * f;
    double n = earthRadius / sqrt(1 - e2 * sin(lat) * sin(lat));
    double x = (n + alt) * cos(lat) * cos(lon);
    double y = (n + alt) * cos(lat) * sin(lon);
    double z = ((1 - e2) * n + alt) * sin(lat);

    Location loc;
    loc.x = x;
    loc.y = y;
    loc.z = z;

    return loc;
}


double SciAtt::computeAngle(Location satellite, Location target1, Position target2) {
        double pi = M_PI;
        double ke1 = satellite.x - target1.x;
        double kn1 = satellite.y - target1.y;
        double ku1 = satellite.z - target1.z;
        double ke2 = ke1 / sqrt(pow(ke1, 2) + pow(kn1, 2) + pow(ku1, 2));
        double kn2 = kn1 / sqrt(pow(ke1, 2) + pow(kn1, 2) + pow(ku1, 2));
        double ku2 = ku1 / sqrt(pow(ke1, 2) + pow(kn1, 2) + pow(ku1, 2));
        double ke = -sin(target2.lon * pi / 180) * ke2 + cos(target2.lon * pi / 180) * kn2 + 0 * ku2;
        double kn = -cos(target2.lon * pi / 180) * sin(target2.lat * pi / 180) * ke2 - sin(target2.lon * pi / 180) * sin(target2.lat * pi / 180) * kn2 + cos(target2.lat * pi / 180) * ku2;
        double ku = cos(target2.lon * pi / 180) * cos(target2.lat * pi / 180) * ke2 + sin(target2.lon * pi / 180) * cos(target2.lat * pi / 180) * kn2 + sin(target2.lat * pi / 180) * ku2;
        double ele = atan(ku / sqrt(pow(ke, 2) + pow(kn, 2))) * 180 / pi;
        double azi = 999.0;
        if (ke > 0 && kn > 0) {
            azi = atan(ke / kn) * 180 / pi;
        } else if (ke >= 0 && kn == 0) {
            azi = 90.0;
        } else if (ke > 0 && kn < 0) {
            azi = atan(ke / kn) * 180 / pi + 180;
        } else if (ke < 0 && kn > 0) {
            azi = atan(ke / kn) * 180 / pi + 360;
        } else if (ke <= 0 && kn == 0) {
            azi = 270;
        } else if (ke < 0 && kn < 0) {
            azi = atan(ke / kn) * 180 / pi + 180;
        }
        return ele;
 }




