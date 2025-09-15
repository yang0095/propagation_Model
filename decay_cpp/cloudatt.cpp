#include "cloudatt.h"


//const double EARTH_RADIUS = 6371.0; // 地球半径 (km)
const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;
//根据经度，纬度查找文件，返回对应Lred Annual Maps文件中的数据
double CloudAtt::findLred(const std::string& longitudeFile, const std::string& latitudeFile, const std::string& LredFile, double inputLongitude, double inputLatitude) {
    std::ifstream lonFile(longitudeFile);
    std::ifstream latFile(latitudeFile);
    std::ifstream lredFile(LredFile);

    if (!lonFile || !latFile || !lredFile) {
        std::cerr << "无法打开其中一个文件。" << std::endl;
        return -1; // 返回-1表示错误
    }

    std::string lonLine, latLine, lredLine;

    // 显式声明lambda表达式的类型为std::function
    std::function<bool(const pair<double, double>&, const pair<double, double>&)> cmp =
        [](const pair<double, double>& p1, const pair<double, double>& p2) {
            return p1.first < p2.first; // 大根堆：第一个元素大的优先
        };

    // 创建一个大根堆，存储pair<double, double>，使用显式声明的lambda表达式,设置大小为2
    priority_queue<pair<double, double>, vector<pair<double, double>>, std::function<bool(const pair<double, double>&, const pair<double, double>&)>> maxHeap(cmp);

    //逐行读取数据，将每组数据和输入的经度纬度求曼哈顿距离，保持最小的两个在大根堆中
    while (std::getline(lonFile, lonLine) && std::getline(latFile, latLine) && std::getline(lredFile, lredLine)) {

        std::istringstream lonStream(lonLine);
        std::istringstream latStream(latLine);
        std::istringstream lredStream(lredLine);

        double longitude, latitude, lred;

        // 逐个读取行内的多个值
        while (lonStream >> longitude && latStream >> latitude && lredStream >> lred) {
            if(longitude==inputLongitude && latitude==inputLatitude){ //如果精确相等，直接返回对应的lred值
                return lred;
            }
            double distance=abs(longitude-inputLongitude)+abs(latitude-inputLatitude);
            pair<double,double>pir=make_pair(distance,lred);
            if(maxHeap.size()<2){
                maxHeap.push(pir);
            }else if(distance<maxHeap.top().first){
                maxHeap.pop();
                maxHeap.push(pir);
            }
        }
    }
    double arv_lred=0.0;
    while(!maxHeap.empty()){
        arv_lred+=maxHeap.top().second;
        maxHeap.pop();
    }
    arv_lred/=2;
    return arv_lred;
}

// 解析文件名中的概率数字
double CloudAtt::parseProbability(const std::string& filename) {
    size_t start = filename.find('_');
    if (start == std::string::npos){
        return -1.0;// 解析失败
    }
    start += 1;

    size_t end = filename.find('_', start);
    if (end == std::string::npos){
        return -1.0;// 解析失败
    }
    std::string numberStr = filename.substr(start, end - start);
    try {
        if(numberStr.size()>1&&numberStr[0]=='0'){
            return std::stod(numberStr)*0.1;
        }else{
            return std::stod(numberStr);
        }
    } catch (const std::invalid_argument&) {
        return -1.0; // 转换失败
    }
}

// 计算与目标概率的差异
double CloudAtt::calculateDifference(double prob1, double prob2) {
    return std::fabs(prob1 - prob2);
}

//根据概率，查找目录下所有Lred文件，选择相等或者最接近的两个文件
std::vector<std::string> CloudAtt::getFileName(const std::string& directoryPath, double targetProbability) {
    std::vector<std::string> result;  // 存放返回的文件名字符串
    std::vector<std::pair<std::string, double>> files; // 存储文件名和对应的概率值

#ifdef _WIN32
    // Windows 下的目录遍历
    std::string searchPath = directoryPath + "\\*";
    WIN32_FIND_DATAA findFileData;
    HANDLE hFind = FindFirstFileA(searchPath.c_str(), &findFileData);

    if (hFind == INVALID_HANDLE_VALUE) {
        std::cerr << "无法打开目录: " << directoryPath << std::endl;
        return result;
    }

    do {
        if (!(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            std::string filename = findFileData.cFileName;
            double probability = parseProbability(filename) * 100;
            if (probability != -1) {
                files.emplace_back(filename, probability);
            }
        }
    } while (FindNextFileA(hFind, &findFileData) != 0);

    FindClose(hFind);
#else
    // Linux 下的目录遍历
    DIR* dir = opendir(directoryPath.c_str());
    if (dir == nullptr) {
        std::cerr << "无法打开目录: " << directoryPath << std::endl;
        return result;
    }
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        if (!(entry->d_type & DT_DIR)) {
            std::string filename = entry->d_name;
            double probability = parseProbability(filename) * 100;
            if (probability != -1) {
                files.emplace_back(filename, probability);
            }
        }
    }

    closedir(dir);
#endif

    // 查找与目标概率最接近的文件
    std::pair<std::string, double> closestFile;
    std::pair<std::string, double> secondClosestFile;
    double minDiff = std::numeric_limits<double>::max();
    double secondMinDiff = std::numeric_limits<double>::max();

    for (const auto& file : files) {
        double diff = calculateDifference(file.second, targetProbability);
        if (diff == 0) {
            result.push_back(file.first);
            return result; // 找到完全相等的概率
        }
        if (diff < minDiff) {
            secondMinDiff = minDiff;
            secondClosestFile = closestFile;
            minDiff = diff;
            closestFile = file;
        } else if (diff < secondMinDiff) {
            secondMinDiff = diff;
            secondClosestFile = file;
        }
    }

    // 输出最接近的文件名
    if (minDiff != std::numeric_limits<double>::max()) {
        result.push_back(closestFile.first);
    } else {
        std::cout << "没有找到最接近的文件。\n";
    }

    // 输出第二接近的文件名
    if (secondMinDiff != std::numeric_limits<double>::max()) {
        result.push_back(secondClosestFile.first);
    } else {
        std::cout << "没有找到第二接近的文件。\n";
    }

    return result;
}


/**
 * 计算云衰减
 * @param f 频率f（Hz） 精度1
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
 * @param cloud_type 气象条件（云）类型
 * @param L 液态水含量(mm) 精度0.10.1
 * @param water_tmp 云中液态水温度（K）精度0.1
 * @return 云衰减（dB）
 */
double CloudAtt::cloud_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                  Receive_polar receive_polar,double dia,double antenna_effeciency, Cloud_type cloud_type,double L,double water_tmp){
    f=f/pow(10,9);//hz->GHz
    if(parameter_rational(f, PP, satlon,  satlat, satalt,  send_polar, lon,  lat, alt, receive_polar, dia, antenna_effeciency,  cloud_type, L, water_tmp)==false){
        exit(0);
    }
    double result = 0.0;
    double t1 = water_tmp;
    double rt = 300.0 / t1;
    double e0 = 77.6 + 103.3 * (rt - 1.0);
    double e1 = 5.48, e2 = 3.51;
    double fp = 20.09 - 142.0 * (rt - 1.0) + 294.0 * pow(rt - 1.0, 2.0);
    double fs = 590.0 - 1500.0 * (rt - 1.0);
    double eppf = f * (e0 - e1) / (fp * (1.0 + pow(f / fp, 2.0))) + f * (e1 - e2) / (fs * (1.0 + pow(f / fs, 2.0)));
    double epf = (e0 - e1) / (1.0 + pow(f / fp, 2.0)) + (e1 - e2) / (1 + pow(f / fs, 2.0)) + e2;
    double yita = (2.0 + epf) / eppf;
    double K = (0.819 * f) / (eppf * (1.0 + pow(yita, 2.0)));

    Position station;
    station.setLat(lat);
    station.setLon(lon);
    station.setAlt(alt);

    Position satellite;
    satellite.setLat(satlat);
    satellite.setLon(satlon);
    satellite.setAlt(satalt);

    Location stationLoc = blhToXyz(station);
    Location satelliteLoc = blhToXyz(satellite);

    double theta = computeAngle(satelliteLoc, stationLoc, station);
    //std::cout<<"仰角："<<theta<<endl;

    result = L * K / sin(theta*M_PI/180);
    double hc=6000.0;//定义云高
    if(satalt>hc){
        if(alt>hc){
            result=0;
        }else{
            double hmin=min(satalt,alt);
            result=result*abs(hc-hmin)/hc;
        }
    }else{
        if(alt>hc){
             double hmin=min(satalt,alt);
            result=result*abs(hc-hmin)/hc;
        }else{
            result=result*abs(alt-satalt)/hc;
        }
    }

    result*=2;
    return std::round(result * 10) / 10;
}

double CloudAtt::specificatten(double f,double L,double water_tmp)
{
//    f=f/pow(10,9);//hz->GHz
    double result = 0.0;
    double t1 = water_tmp;
    double rt = 300.0 / t1;
    double e0 = 77.6 + 103.3 * (rt - 1.0);
    double e1 = 5.48, e2 = 3.51;
    double fp = 20.09 - 142.0 * (rt - 1.0) + 294.0 * pow(rt - 1.0, 2.0);
    double fs = 590.0 - 1500.0 * (rt - 1.0);
    double eppf = f * (e0 - e1) / (fp * (1.0 + pow(f / fp, 2.0))) + f * (e1 - e2) / (fs * (1.0 + pow(f / fs, 2.0)));
    double epf = (e0 - e1) / (1.0 + pow(f / fp, 2.0)) + (e1 - e2) / (1 + pow(f / fs, 2.0)) + e2;
    double yita = (2.0 + epf) / eppf;
    double K = (0.819 * f) / (eppf * (1.0 + pow(yita, 2.0)));
    //(zr)确认？
    result=K*L;
    return result;
}

double CloudAtt::FUN_cloud_att(double f,double lat1, double lon1, double h1,
                               double lat2, double lon2, double h2,
                               int cloud_area_count,vector<vector<double>>cloud_area,vector<double>cR)
{

    double water_tmp=273.15;
    double cloud_att_count=0.0;
    for(int i=0;i<cloud_area_count;i++){
        if(parameter_rational1(f,cR[i],water_tmp,lat1, lon1,lat2, lon2,
                               cloud_area[i][0], cloud_area[i][1], cloud_area[i][3], cloud_area[i][4])==false){
            exit(0);
        }
        try {
                // 计算有效路径长度
                double r = calculatePathLength(
                    lat1, lon1, h1, lat2, lon2, h2,
                    cloud_area[i][0], cloud_area[i][1], cloud_area[i][3], cloud_area[i][4],
                    cloud_area[i][2], cloud_area[i][5]);

                std::cout << "云区内有效路径长度: " << r << " 米\n";

                // 计算衰减
                double atten_km=specificatten(f,  cR[i], water_tmp);
                std::cout << "云衰减率: " << atten_km ;
                cloud_att_count+= atten_km * r/1000;

            } catch (const std::exception& e) {
                std::cerr << "错误: " << e.what() << std::endl;
                return 1;
            }
    }
    return cloud_att_count;

}

double CloudAtt::calculatePathLength(double lat1, double lon1, double h1, double lat2, double lon2, double h2, double lata, double lona, double latb, double lonb, double ha, double hb)
{
    CcalculateElevation mycalelev;
    cisPlaneBetween myplanbtn;
    cisPlaneBetween m_plantn;
    CdegTorad mydegrad;
    CedgeRectangle1 myedgrect1;
    CedgeRectangle2 myedgrect2;
    CcalculateHeight mycalheight;
    // 确保雷达高度 < 目标高度
    cout<<"初始H1 H2 Ha Hb"<<h1<<" "<<h2<<" "<<ha<<" "<<hb<<endl;
       if (h1 > h2) {
           std::swap(h1, h2);
       }

       // 确保云底高度 < 云顶高度
       if (ha > hb) {
           std::swap(ha, hb);
       }

       // 计算仰角
       double elev_deg = mycalelev.calculateElevation(lat1,lon1,h1,lat2,lon2,h2);
       cout<<"计算得到角度(云)"<<elev_deg<<endl;
       double elev_rad = elev_deg * DEG_TO_RAD;
       // 高度分类 (ih)
       int ih = 0;
       if (h2 < ha) ih = 1;
       else if (h1 < ha) {
           ih = (h2 < hb) ? 2 : 3;
       }
       else if (h2 < hb) ih = 4;
       else if (h1 < hb) ih = 5;
       else ih = 6;

       // 平面位置分类 (ip)
       int ip = 0;
       bool b1,b2;
       if(lat1==lata ||lat1==latb || lon1==lona ||lon1==lonb  ){ b1 = myplanbtn.isPlaneBetween(lat1,lon1,lata,lona,latb,lonb,true);}
       else {
           b1 = myplanbtn.isPlaneBetween(lat1,lon1,lata,lona,latb,lonb);
       }
       if(lat2==lata || lat2==latb|| lon2==lona || lon2==lonb){ b2 = myplanbtn.isPlaneBetween(lat2,lon2,lata,lona,latb,lonb,true);}
       else {
            b2 = myplanbtn.isPlaneBetween(lat2,lon2,lata,lona,latb,lonb);
       }

       if (b1 && b2) ip = 1;
       else if (b1) ip = 2;
       else if (b2) ip = 3;
       else ip = 4;

       cout<<"B1 & B2 ::"<<b1<<" "<<b2<<endl;

       cout<<"Ip & Ih"<<ip<<" "<<ih<<endl;

       // 计算有效路径长度 r
       double r = 0.0;
       double dh = 0.0;
       switch (ih) {
           case 2: {//地面基站高度小于云底，飞行器高度大于云底小于云顶
           cout<<"ih is 2"<<endl;
               switch (ip) {
                   case 1:{
                   dh = h2 -ha;
                   r = dh / sin(elev_rad);
                   break;
               }

               case 2:{
                   CedgeRectangle2 myrectangle2;
                   CcalculateHeight mycalheight;
                   auto [lat, lon] = myrectangle2.edgeRectangle2(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
                  cout<<lat.size()<<" "<<lat[0]<<endl;
                   if (lat.size() == 2) {
                       double h1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                       double h2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[1], lon[1]);
                       dh = std::abs(h1 - h2);
                       h1 = std::max(h1, ha);
                       h2 = std::max(h2, ha);
                       r = dh / sin(elev_rad);
                   }
               break;
               }
                   case 3: {

                       CedgeRectangle1 myrectangle1;
                       CcalculateHeight mycalheight;
                       auto [lat, lon] = myrectangle1.edgeRectangle1(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
                       if (!lat.empty()) {
                           double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                           hcrs = std::max(hcrs, ha);
                           dh = (ip == 2) ? (hcrs - ha) : (h2 - hcrs);
                           r = dh / sin(elev_rad);
                       }
                       break;
                   }
                   case 4: {

                       CedgeRectangle2 myrectangle2;
                       CcalculateHeight mycalheight;
                       auto [lat, lon] = myrectangle2.edgeRectangle2(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
                      cout<<lat.size()<<endl;
                       if (lat.size() == 2) {
                           double h1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                           double h2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[1], lon[1]);
                           dh = std::abs(h1 - h2);
                           h1 = std::max(h1, ha);
                           h2 = std::max(h2, ha);
                           r = dh / sin(elev_rad);

//                           if (h1 >= ha || h2 >= ha) {
//                               h1 = std::max(h1, ha);
//                               h2 = std::max(h2, ha);
//                               dh = std::abs(h1 - h2);
//                               r = dh / sin(elev_rad);
//                           }
                       }
                       break;
                   }
               }
               break;
           }
      case 3:{//地面雷达高度小于云底，飞行器高度大于云顶
           cout<<"ih is 3"<<endl;
           switch (ip) {
               case 1: // 两点都在矩形内
           {
               r = (hb - ha) / sin(elev_rad);
               cout<<"h2_h1"<<h2<<" "<<h1<<sin(elev_rad)<<endl;
               break;
           }


               case 2: { // 站点1在矩形内，站点2在外
               CedgeRectangle2 myrectangle2;
                   auto [lat, lon] = myedgrect1.edgeRectangle1(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
                   if (!lat.empty()) {
                       double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                       hcrs = std::min(ha, hcrs);
                       cout<<hcrs<<"is hcrs"<<endl;
                       cout<<elev_deg<<" "<<sin(elev_rad)<<endl;
                       r = (hcrs - h1) / sin(elev_rad);
                   }
                   break;
               }

           case 3: {

               CedgeRectangle1 myrectangle1;
               CcalculateHeight mycalheight;
               auto [lat, lon] = myrectangle1.edgeRectangle1(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
               if (!lat.empty()) {
                   double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                   hcrs = std::max(hcrs, ha);
                   dh = (ip == 2) ? (hcrs - ha) : (h2 - hcrs);
                   r = dh / sin(elev_rad);
               }
               break;
           }

           case 4: { // 两点都在矩形外
               CedgeRectangle2 myrectangle2;
               CcalculateHeight mycalheight;
               auto [lat, lon] = myrectangle2.edgeRectangle2(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
              cout<<lat.size()<<endl;
               if (lat.size() == 2) {
                   double h1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                   double h2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[1], lon[1]);
                   dh = std::abs(h1 - h2);
                   h1 = std::max(h1, ha);
                   h2 = std::max(h2, ha);
                   r = dh / sin(elev_rad);
               }
               break;
           }
       }
           break;
    }

       case 4:{//地面基站高度大于云底，飞行器高度大于云底小于云顶
           cout<<"ih is 4"<<endl;
           switch (ip){
           case 1:{//两点都在矩形内，直接计算
               r = (h2 - h1) / sin(elev_rad);
               break;
           }
           case 2: { // 站点1在矩形内，站点2在外
               auto [lat, lon] = myedgrect1.edgeRectangle1(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
               if (!lat.empty()) {
                   double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                   r = (hcrs - h1) / sin(elev_rad);
                   cout<<"hcrs h1 r"<<hcrs<<" "<<h1<<" "<<r<<" "<<sin(elev_rad)<<endl;
                   break;
               }

               break;
           }
           case 3: { // 站点2在矩形内，站点1在外
               auto [lat, lon] = myedgrect1.edgeRectangle1(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
               if (!lat.empty()) {
                   double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
//                        cout<<"计算高度为"<<hcrs<<endl;
                   r = (h2 - hcrs) / sin(elev_rad);
               }
               break;
           }

           case 4: { // 两点都在矩形外
               auto [latn, lonn] = myedgrect2.edgeRectangle2(lat1, lon1, lat2, lon2, lata, lona, latb, lonb);
//                            cout<<"latn size is"<<latn.size()<<endl;
               if (latn.size() == 2) {
                   double hcrs1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, latn[0], lonn[0]);
                   double hcrs2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, latn[1], lonn[1]);
                   r = std::abs(hcrs1 - hcrs2) / sin(elev_rad);
               }
               break;
           }

           }
          }
       case 5: {//地面基站高度大于云底，飞行器高度大于云顶
           switch (ip) {//两点都在范围内
               case 1:
                   dh = hb - h1;
                   r = dh / sin(elev_rad);
                   break;

           case 2:{
                   CedgeRectangle2 myrectangle2;
                   CcalculateHeight mycalheight;
                  auto [lat, lon] = myrectangle2.edgeRectangle2(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
                  if (!lat.empty()) {
                      double hcrs1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                      double hcrs2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[1], lon[1]);
                      dh = abs(hcrs1-hcrs2);
                      r = dh / sin(elev_rad);
                  }
                 break;

           }
           case 3: {

                   CedgeRectangle1 myrectangle1;
                   CcalculateHeight mycalheight;
                   auto [lat, lon] = myrectangle1.edgeRectangle1(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
                   if (!lat.empty()) {
                       double hcrs = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                       hcrs = std::max(hcrs, ha);
                       dh = (ip == 2) ? (hcrs - ha) : (h2 - hcrs);
                       r = dh / sin(elev_rad);
                   }
                   break;
               }
               case 4: {

                   CedgeRectangle2 myrectangle2;
                   CcalculateHeight mycalheight;
                   auto [lat, lon] = myrectangle2.edgeRectangle2(lat1,lon1,lat2,lon2,lata,lona,latb,lonb);
                  cout<<lat.size()<<endl;
                   if (lat.size() == 2) {
                       double h1 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[0], lon[0]);
                       double h2 = mycalheight.calculateHeight(lat1, lon1, h1, lat2, lon2, h2, lat[1], lon[1]);
                       dh = std::abs(h1 - h2);
                       h1 = std::max(h1, ha);
                       h2 = std::max(h2, ha);
                       r = dh / sin(elev_rad);

//                           if (h1 >= ha || h2 >= ha) {
//                               h1 = std::max(h1, ha);
//                               h2 = std::max(h2, ha);
//                               dh = std::abs(h1 - h2);
//                               r = dh / sin(elev_rad);
//                           }
                   }
                   break;
               }
           }
           break;
       }



           // 其他情况类似处理 (ih=3,4,5)
           // ...

           default:
               // 处理其他情况
               break;
       }

       // 输出结果
       std::cout << "Effective path length: " << r << " meters" << std::endl;
       return r;
}


bool CloudAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                        Receive_polar receive_polar,double dia,double antenna_effeciency, Cloud_type cloud_type,double L,double water_tmp)
{

    bool flag=true;
    if(f<12||f>18){
        std::cout<<"cloudatt:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"cloudatt:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"cloudatt:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"cloudatt:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"cloudatt:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"cloudatt:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(dia<0||dia>5){
        std::cout<<"天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
        exit(0);
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }
    if(L<0||L>10){
        std::cout<<"cloudatt:液态水含量="<<L<<"超出限定范围（0~10）!"<<endl;
        flag=false;
    }

    if(water_tmp<193.5||water_tmp>323.5){
        std::cout<<"cloudatt:云中液态水温度="<<water_tmp<<"超出限定范围（193.5~323.5）!"<<endl;
         flag=false;
    }
    return flag;
}

bool CloudAtt::parameter_rational1(double f, double L, double water_tmp, double lat1, double lon1, double lat2, double lon2, double lata, double lona, double latb, double lonb)
{
    bool flag=true;
    if(f<12||f>18){
        std::cout<<"cloudatt:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }
    if(L<0||L>10){
        std::cout<<"cloudatt:液态水含量="<<L<<"超出限定范围（0~10）!"<<endl;
        flag=false;
    }

    if(water_tmp<193.5||water_tmp>323.5){
        std::cout<<"cloudatt:云中液态水温度="<<water_tmp<<"超出限定范围（193.5~323.5）!"<<endl;
         flag=false;
    }
    if(lon1<-180||lon1>180){
        std::cout<<"cloudatt:雷达经度="<<lon1<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat1<-90||lat1>90){
        std::cout<<"cloudatt:雷达纬度="<<lat1<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(lon2<-180||lon2>180){
        std::cout<<"cloudatt:目标经度="<<lon2<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat2<-90||lat2>90){
        std::cout<<"cloudatt:目标纬度="<<lat2<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(lona<-180||lona>180){
        std::cout<<"cloudatt:左下经度="<<lona<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lata<-90||lata>90){
        std::cout<<"cloudatt:左下纬度="<<lata<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(lonb<-180||lonb>180){
        std::cout<<"cloudatt:右上经度="<<lonb<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(latb<-90||latb>90){
        std::cout<<"cloudatt:右上纬度="<<latb<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    return flag;
}

Location CloudAtt::blhToXyz(Position p) {
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


double CloudAtt::computeAngle(Location satellite, Location target1, Position target2) {
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



