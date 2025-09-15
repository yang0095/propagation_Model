#include "rainatt.h"



//根据经度，维度查找文件，返回对应降雨率文件中的数据
double RainAtt::findValue(const std::string& longitudeFile, const std::string& latitudeFile, const std::string& nWetFile, double inputLongitude, double inputLatitude) {
    std::ifstream lonFile(longitudeFile);
    std::ifstream latFile(latitudeFile);
    std::ifstream wetFile(nWetFile);

    if (!lonFile || !latFile || !wetFile) {
        std::cerr << "无法打开其中一个文件。" << std::endl;
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
    return std::round(arv_n_wet * 10) / 10;
}


/**
 * @brief CalAttRain
 * @param f 频率f（Hz->GHz） 精度1
 * @param percentage 时间概率（%） 精度0.001
 * @param satlon 卫星经度（度）精度10-6
 * @param satlat 卫星纬度（度）精度10-6
 * @param satalt 飞行器高度（m）精度0.1
 * @param send_polar  信源发射天线极化方式
 * @param site_lon 站点经度（度）精度10-6
 * @param site_lat 站点纬度（度）精度10-6
 * @param altitude 站点海拔高度h_s（m）精度0.1
 * @param receive_polar  信宿接收天线极化方式
 * @param dia 直径（m） 精度0.1
 * @param antenna_effeciency 天线效率 精度0.01
 * @param meteor_condition 气象条件等级
 * @param rainrate 降雨率
 * @param start_distance 气象条件（雨）降雨起始距离m 精度0.1
 * @param end_distance 气象条件（雨）降雨终止距离m 精度0.1
 * @param rain_height 雨顶高度m 精度0.1
 */
double RainAtt::rain_att(double f,double percentage,double satlon, double satlat,double satalt, Send_polar send_polar,double site_lon, double site_lat,double altitude,
                 Receive_polar receive_polar,double dia,double antenna_effeciency, Meteor_condition meteor_condition,double rainrate,
                double start_distance,double end_distance,double rain_height){
    //判断参数合理性
    if(parameter_rational( f, percentage, satlon,  satlat, satalt,  send_polar, site_lon,  site_lat, altitude,receive_polar, dia, antenna_effeciency, rainrate, rain_height)==false){
        exit(0);
    }

    // 计算雨衰减
     f=f/pow(10,9);//HZ->GHz
    Position station;
    station.setLat(site_lat);
    station.setLon(site_lon);
    station.setAlt(altitude);

    Position satellite;
    satellite.setLat(satlat);
    satellite.setLon(satlon);
    satellite.setAlt(satalt);

    Location stationLoc = blhToXyz(station);
    Location satelliteLoc = blhToXyz(satellite);

    double elevation = computeAngle(satelliteLoc, stationLoc, station);

    int tmp=static_cast<int>(meteor_condition);
    if(0==tmp){
        rainrate=std::max(rainrate,4.0);
    }else if(1==tmp){
         rainrate=std::max(rainrate,14.0);
    }else{
         rainrate=std::max(rainrate,40.0);
    }

    double tao;
    if(static_cast<int>(send_polar)==0||static_cast<int>(send_polar)==1){
        tao=0;
    }else{
        tao=90;
    }

    double frequency=f;
    std::vector<double>temp=specificatten(elevation, tao, frequency, rainrate);
    double atten_km=temp[2];
    double alpha=temp[1];

    double RainHeight =std::min(rain_height,max(altitude,satalt))/1000  + 0.36;
    double SlantLength = GetSLength(elevation, altitude/1000, RainHeight);
    double GroundLength = SlantLength * cos(elevation * M_PI / 180.0);
    double xlh = 0.78 * sqrt(GroundLength * atten_km / frequency);
    double qmx = 0.38 * (1 - exp(-2 * GroundLength));
    double HoriFactor = 1 / (1 + xlh - qmx);
    double cita = atan((RainHeight - altitude) / (GroundLength * HoriFactor)) * 180 / M_PI;
    double LR;
    if (cita > elevation) {
        LR = GroundLength * HoriFactor / cos(elevation * M_PI / 180.0);
    } else {
        LR = (RainHeight - altitude/1000) / sin(elevation * M_PI / 180.0);
    }

    double x;
    if (fabs(site_lat) < 36) {
        x = 36 - fabs(site_lat);
    } else {
        x = 0;
    }
    double lcs = sqrt(sin(elevation * M_PI / 180.0));
    double lx = (1 - exp(-(elevation / (1 + x)))) * sqrt(LR * atten_km) / (frequency * frequency);
    double v01 = 1 / (1 + lcs * (31 * lx - 0.45));
    double EffectiveLength = LR * v01;
    double PredictAtten01 = atten_km * EffectiveLength;

    double prediction;

    if (percentage == 0.01) {
        prediction = PredictAtten01;
    } else {
        double beita, coeff_mie;
        if (percentage >= 1 && fabs(site_lat) >= 25) {
            beita = 0;
        } else if (percentage < 1 || elevation >= 25 || fabs(site_lat) >= 36) {
            beita = -0.005 * (fabs(site_lat) - 36);
        } else {
            beita = -0.005 * (fabs(site_lat) - 36) + 1.8 - 4.25 * sin(cita * M_PI / 180.0);
        }
        coeff_mie = -(0.655 + 0.033 * log(percentage)) - 0.045 * log(PredictAtten01) - beita * (1 - percentage) * sin(cita * M_PI / 180.0);
        prediction = PredictAtten01 * pow(percentage / 0.01, coeff_mie);
    }

    //（zr）拷贝过去
    if(1==static_cast<int>(send_polar)||2==static_cast<int>(send_polar)){
        double xpd_value=xpd(frequency,percentage,elevation,tao,prediction);
        prediction+=prediction+xpd_value-3;
    }else{
        prediction*=2;
    }
   return std::round(prediction * 10) / 10;
}

Location RainAtt::blhToXyz(Position p) {
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


double RainAtt::computeAngle(Location satellite, Location target1, Position target2) {
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


 bool RainAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                  Receive_polar receive_polar,double dia,double antenna_effeciency,double rainrate,double rain_height){

     bool flag=true;
     f=f/pow(10,9);
     if(f<12||f>18){
         std::cout<<"rain_att:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
         flag=false;
     }

     if(satlon<-180||satlon>180){
         std::cout<<"rain_att:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
         flag=false;
     }
     if(satlat<-90||satlat>90){
         std::cout<<"rain_att:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
         flag=false;
     }
     double h=satlat/1000;
     if(h<0||h>36000){
         std::cout<<"rain_att:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
         flag=false;
     }
     if(lon<-180||lon>180){
         std::cout<<"rain_att:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
         flag=false;
     }
     if(lat<-90||lat>90){
         std::cout<<"rain_att:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
         flag=false;
     }
     if(dia<0||dia>5){
         std::cout<<"rain_att:天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
         exit(0);
     }
     if(antenna_effeciency<0||antenna_effeciency>1){
         std::cout<<"rain_att:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
         exit(0);
     }

     if(rainrate<0||rainrate>150){
         std::cout<<"rain_att:降雨率="<<rainrate<<" 不符合范围（0~150)"<<std::endl;
         exit(0);
     }

     if(rain_height<0||rain_height>20000){
         std::cout<<"rain_att:雨顶高度="<<rain_height<<" 不符合范围（0~20000)m"<<std::endl;
         exit(0);
     }
     return flag;
 }

 double RainAtt::calculateRainPathLength(double lat1, double lon1, double h1, double lat2, double lon2, double h2, double lata, double lona, double latb, double lonb, double hR, double elev_d)
 {
     cisPlaneBetween m_plantn;
     CdegTorad mydegrad;
     CedgeRectangle1 myedgrect1;
     CedgeRectangle2 myedgrect2;
     CcalculateHeight mycalheight;
     // 高度分类
        int ih = 0;
        if (h2 < hR) {
            ih = 1;
        } else if (h1 < hR) {
            ih = 2;
        } else {
            ih = 3;
        }

        // 平面分类
        int ip = 0;
        bool b1 = m_plantn.isPlaneBetween(lat1, lon1, lata, lona, latb, lonb);
        bool b2 = m_plantn.isPlaneBetween(lat2, lon2, lata, lona, latb, lonb);


        if (b1 && b2) {
            ip = 1;
        } else if (b1) {
            ip = 2;
        } else if (b2) {
            ip = 3;
        } else {
            ip = 4;
        }

        cout<<"Ip & Ih"<<ip<<" "<<ih<<endl;

        // 计算仰角弧度
        double elev_r = mydegrad.deg2rad(elev_d);
        double r = 0.0; // 有效路径长度
//        cout<<"IH is:"<<ih<<"\n"<<"IP is"<<ip<<"\n"<<"b1 & b2 is"<<b1<<" "<<b2<<endl;

        // 分情况计算
        if (ih == 1) {
            switch (ip) {
                case 1: // 两点都在矩形内    NEW:直接计算距离
                    r = (h2 - h1) / sin(elev_r);
                    cout<<"h2 is "<<h2<<"h1 is "<<h1<<"ele is"<<sin(-elev_r)<<"R is"<<r<<endl;

                    break;

                case 2: { // 站点1在矩形内，站点2在外
                    auto [lat, lon] = myedgrect1.edgeRectangle1(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
                    if (!lat.empty()) {
                        double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                        r = (hcrs - h1) / sin(elev_r);
                        cout<<"hcrs h1 r"<<hcrs<<" "<<h1<<" "<<r<<endl;
                    }
                    break;
                }

                case 3: { // 站点2在矩形内，站点1在外
                    auto [lat, lon] = myedgrect1.edgeRectangle1(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
                    if (!lat.empty()) {
                        double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
//                        cout<<"计算高度为"<<hcrs<<endl;
                        r = (h2 - hcrs) / sin(elev_r);
                    }
                    break;
                }

                case 4: { // 两点都在矩形外
                    auto [latn, lonn] = myedgrect2.edgeRectangle2(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
//                            cout<<"latn size is"<<latn.size()<<endl;
                    if (latn.size() == 2) {
                        double hcrs1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, latn[0], lonn[0]);
                        double hcrs2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, latn[1], lonn[1]);
                        r = std::abs(hcrs1 - hcrs2) / sin(elev_r);
                    }
                    break;
                }
            }
        } else if (ih == 2) {
            switch (ip) {
                case 1: // 两点都在矩形内
                    r = (hR - h1) / sin(elev_r);
                    break;

                case 2: { // 站点1在矩形内，站点2在外
                    auto [lat, lon] = myedgrect1.edgeRectangle1(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
                    if (!lat.empty()) {
                        double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                        hcrs = std::min(hR, hcrs);
                        r = (hcrs - h1) / sin(elev_r);
                    }
                    break;
                }

                case 3: { // 站点2在矩形内，站点1在外
                    auto [lat, lon] = myedgrect1.edgeRectangle1(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
                    if (!lat.empty()) {
                        double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                        hcrs = std::min(hR, hcrs);
                        //cout<<"计算高度为"<<hcrs<<endl;
                        r = (hR - hcrs) / sin(elev_r);
                    }
                    break;
                }

                case 4: { // 两点都在矩形外
                    auto [latn, lonn] = myedgrect2.edgeRectangle2(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
                    if (latn.size() == 2) {
                        double hcrs1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, latn[0], lonn[0]);
                        double hcrs2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, latn[1], lonn[1]);
                        hcrs1 = std::min(hR, hcrs1);
                        hcrs2 = std::min(hR, hcrs2);
                        r = std::abs(hcrs1 - hcrs2) / sin(elev_r);
                    }
                    break;
                }
            }
        }

        // 确保路径长度非负
        if(r<=0) r=-r;
        return std::max(0.0, r);
 }
//zy 这里添加一个mX6的矩阵，代表m组雨区的左下角经纬高、右上角经纬高  vector<vecotr<6>>
//zy 额外添加一个mX1的矩阵，代表m组雨区的降雨类型    elev_d:仰角 rR：降雨率  frequency：频率   polarization：极化方式
 double RainAtt::FUN_rain_att(double lat1, double lon1, double h1, double lat2, double lon2, double h2,
                              int rain_area_count,vector<vector<double>>rain_area,
                              vector<double > rR, double frequency, int polarization)
 {
     try {
             double att=0.0;
             double elev_d;


             Position satellite;
             satellite.setLat(lat1);
             satellite.setLon(lon1);
             satellite.setAlt(h1);

             Position station;
             station.setLat(lat2);
             station.setLon(lon2);
             station.setAlt(h2);

             Location satelliteLoc = blhToXyz(satellite);
             Location stationLoc = blhToXyz(station);

             elev_d = computeAngle(satelliteLoc, stationLoc, station);
             cout<<"计算角度为"<<elev_d<<endl;



             // 计算有效路径长度
             for(int i=0;i<rain_area_count;i++){
 //                cout<<"HR::"<<rain_area[i][5]<<endl;
                 double r = calculateRainPathLength (
                     lat1, lon1, h1, lat2, lon2, h2,
                     rain_area[i][0], rain_area[i][1],
                     rain_area[i][3], rain_area[i][4],
                     rain_area[i][5], elev_d
                 );

                 std::cout <<"The no. "<<i+1<<" Rain Area effic_Length is : " << r << " m"<<endl;

                 std::vector<double>temp=specificatten(elev_d, polarization, frequency, rR[i]);  //计算衰减

                 std::cout <<"The no. "<<i+1<<" Rain Area att is : " << temp[2] << " "<<endl;

                 att+=r*temp[2]/1000;
             }
             std::cout <<"Total Rain Area  att  : " << att <<endl;
             return att;

         } catch (const std::exception& e) {
             std::cerr << "错误: " << e.what() << std::endl;
             return 1;
         }
         return 0;
 }

 //             double r = calculateRainPathLength (
 //                 lat1, lon1, h1, lat2, lon2, h2,
 //                 lata, lona, latb, lonb,
 //                 hR, elev_d
 //             );



              // 计算雨衰减
 //             std::vector<double>temp=specificatten(elev_d, polarization, frequency, rR);  //(zr)发射极化
              //zy1.计算发射极化
 //             double atten_km=temp[2];
 //             std::cout << "雨衰减率: " << atten_km ;
 //             return atten_km * r;
              //(zr)参数+极化参数；计算交叉极化




//%函数功能：计算降雨的交叉极化鉴别度cross-polarization discrimination
//% 适用条件6<=f<=55, th<=60;
//%输入参数：
//% f——频率,GHz;
//% p——时间概率，%；
//% th——路径仰角，度；
//% tau——极化角，度
//% ap——降雨的同极化衰减，dB.

double RainAtt::xpd(double f,double p,double th,double tau,double ap){
    double cf=0;
    if(f>=6&&f<9){
       cf = 60*log10(f)-28.3;
    }else if(f<36){
        cf = 26*log10(f)+4.1;
    }else if(f<=55){
        cf = 35.9*log10(f)-11.3;
    }

    double vf=0;
    if (f>=6 && f<9){
        vf = 30.8*pow(f,-0.21);
    }else if(f<20){
        vf = 12.8*pow(f,0.19);
    }else if( f<40){
        vf = 22.6;
    }else if(f<55){
        vf = 13.0*pow(f,0.15);
    }

    double ca = vf*log10(ap);
    double ct = -10*log10( 1-0.484*( 1+cos(4*tau*M_PI/180) ) );
    double cth = -40*log10( cos(th*M_PI/180) );

    std::vector<double>pn{1, 0.1, 0.01, 0.01};
    std::vector<double>sign{0, 5, 10, 15};
   double sigma=nearestExtrapolate(pn,sign,p);
   //sigma = interp1(pn,sign,p,'nearest','extrap'); %(zr)先暂时这么用吧
   double csig = 0.0053*sigma*sigma;
   double xpdr = cf - ca + ct +cth+csig;
   double ci = xpdr* (0.3+0.1*log10(p))/2;
   double xpdp = xpdr-ci;
   return xpdp;
}

double RainAtt::nearestExtrapolate(const std::vector<double>&pn,const std::vector<double>&sign,double p){
    if(pn.empty()||sign.empty()||pn.size()!=sign.size()){
        throw std::invalid_argument("输入的数组为空或者长度不相等");
    }
    if(p<pn.front()){
        return sign.front();
    }else if(p>pn.back()){
        return sign.back();
    }

    double closestValue=sign.front();
    double minDistance=std::numeric_limits<double>::max();
    for(size_t i=0;i<pn.size();++i){
        double distance=std::abs(pn[i]-p);
        if(distance<minDistance){
            closestValue=sign[i];
        }
    }
    return closestValue;
}



// 实现雨衰减率的计算
/**
     * @param theta 仰角
     * @param tao   极化方式，圆极化（45），水平极化（0），垂直极化（90）
     * @param f     频率（GHz）
     * @param r     降雨率
     * @return [0]kk [1]aa [2]att雨衰减率 dB/km
     */
std::vector<double> RainAtt::specificatten(double theta, double tao, double f, double r) {
    std::vector<double> result(3);
    double kh = 0, kv = 0, ah = 0, av = 0;
    std::vector<double> a_kh = {-5.3398, -0.35351, -0.23789, -0.94158};
    std::vector<double> b_kh = {-0.10008, 1.2697, 0.86036, 0.64552};
    std::vector<double> c_kh = {1.13098, 0.454, 0.15354, 0.16817};
    double mk_kh = -0.18961;
    double ck_kh = 0.71147;

    std::vector<double> a_ah = {-0.14318, 0.29591, 0.32177, -5.3761, 16.1721};
    std::vector<double> b_ah = {1.82442, 0.77564, 0.63773, -0.9623, -3.2998};
    std::vector<double> c_ah = {-0.55187, 0.19822, 0.13164, 1.47828, 3.4399};
    double ma_ah = 0.67849;
    double ca_ah = -1.95537;

    std::vector<double> a_kv = {-3.80595, -3.44965, -0.39902, 0.50167};
    std::vector<double> b_kv = {0.56934, -0.22911, 0.73042, 1.07319};
    std::vector<double> c_kv = {0.81061, 0.51059, 0.11899, 0.27195};
    double mk_kv = -0.16398;
    double ck_kv = 0.63297;

    std::vector<double> a_av = {-0.07771, 0.56727, -0.20238, -48.2991, 48.5833};
    std::vector<double> b_av = {2.3384, 0.95545, 1.1452, 0.791669, 0.791459};
    std::vector<double> c_av = {-0.76284, 0.54039, 0.26809, 0.116226, 0.116479};
    double ma_av = -0.053739;
    double ca_av = 0.83433;

    double sum_kh = 0, sum_kv = 0, sum_ah = 0, sum_av = 0;

    for (int i = 0; i < 4; ++i) {
        sum_kh += a_kh[i] * std::exp(-std::pow((std::log10(f) - b_kh[i]) / c_kh[i], 2));
        sum_kv += a_kv[i] * std::exp(-std::pow((std::log10(f) - b_kv[i]) / c_kv[i], 2));
    }
    kh = std::pow(10, sum_kh + mk_kh * std::log10(f) + ck_kh);
    kv = std::pow(10, sum_kv + mk_kv * std::log10(f) + ck_kv);

    for (int j = 0; j < 5; ++j) {
        sum_ah += a_ah[j] * std::exp(-std::pow((std::log10(f) - b_ah[j]) / c_ah[j], 2));
        sum_av += a_av[j] * std::exp(-std::pow((std::log10(f) - b_av[j]) / c_av[j], 2));
    }
    ah = sum_ah + ma_ah * std::log10(f) + ca_ah;
    av = sum_av + ma_av * std::log10(f) + ca_av;

    double k = (kh + kv + (kh - kv) * std::cos(theta * M_PI / 180) * std::cos(theta * M_PI / 180) * std::cos(2 * tao * M_PI / 180)) / 2;
    double alpha = (kh * ah + kv * av + (kh * ah - kv * av) * std::cos(theta * M_PI / 180) * std::cos(theta * M_PI / 180) * std::cos(2 * tao * M_PI / 180)) / (2 * k);
    double att = k * std::pow(r, alpha);

    result[0] = k;
    result[1] = alpha;
    result[2] = att;

    return result;
}

// 计算斜路径长度
double RainAtt::GetSLength(double theta, double altitude, double RainHeight) {
    const double RADIU_EARTH = 8500.0; // 等效地球半径
    if (theta < 5) {
        return 2 * (RainHeight - altitude) / (sin(theta * M_PI / 180.0) +
                                              sqrt(pow(sin(theta * M_PI / 180.0), 2) +
                                                   2 * (RainHeight - altitude) / RADIU_EARTH));
    } else {
        return (RainHeight - altitude) / sin(theta * M_PI / 180.0);
    }
}

//zy1 计算发射极化参数
double RainAtt::GetPolar(double polar){
    double tao;
    if(static_cast<int>(polar)==0||static_cast<int>(polar)==1){
        tao=0;
    }else{
        tao=90;
    }

}

