#ifndef CLOUDATT_H
#define CLOUDATT_H
#include "reconfig.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include<iomanip>
#include<ctime>
#include<queue>
#include<functional>
#include <cstring>
#include <string>
#include <utility>
#include <cmath>
#include <limits>
#include <cstdlib>
#include"enum.h"

#include"ccalculateelevation.h"
#include"cisplanebetween.h"
#include"cedgerectangle1.h"
#include"cedgerectangle2.h"
#include"ccalculateheight.h"

#include"rainatt.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <cstring>
#endif
using namespace std;


class  CloudAtt
{  
    public:

    /**
     * 计算云衰减
     * @param f 频率f（单位：Hz） 精度1
     * @param PP 时间概率（单位：%） 精度0.001
     * @param satlon 卫星经度（单位：度）精度10-6
     * @param satlat 卫星纬度（单位：度）精度10-6
     * @param satalt 飞行器高度（单位：m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param lon 站点经度（单位：度）精度10-6
     * @param lat 站点纬度（单位：度）精度10-6
     * @param alt 站点高度（单位：m) 精度0.1
     * @param receive_polar  信宿接收天线极化方式
     * @param dia 直径（单位：m） 精度0.1
     * @param antenna_effeciency 天线效率 精度0.01
     * @param cloud_type 气象条件（云）类型
     * @param L 液态水含量(单位：mm) 精度0.10.1
     * @param water_tmp 云中液态水温度（单位：K）精度0.1
     * @return 云衰减
     */
    double cloud_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                      Receive_polar receive_polar,double dia,double antenna_effeciency, Cloud_type cloud_type,double L,double water_tmp);


    /**
     * 计算云衰减率
     * @param f 频率f（单位：Hz） 精度1
     * @param L 液态水含量(单位：mm) 精度0.10.1
     * @param water_tmp 云中液态水温度（单位：K）精度0.1
     * @return 云衰减率
     */

    double specificatten(double f,double L,double water_tmp);

    /**
     * @brief 计算衰减
     * @param f 频率f（单位：Hz） 精度1
     * @param L 液态水含量(单位：mm) 精度0.10.1
     * @param water_tmp 云中液态水温度（单位：K）精度0.1
     * @param lat1, lon1, h1 雷达的纬度、经度、高度
     * @param lat2, lon2, h2 目标的纬度、经度、高度
     * @param lata, lona, latb, lonb 左下右上两个对角点的纬度和经度
     * @param ha，hb 云底，云顶
     * @return 衰减值（dB）
     */
    double FUN_cloud_att(double f,double lat1, double lon1, double h1,
                         double lat2, double lon2, double h2,
                         int cloud_area_count,vector<vector<double>>cloud_area,vector<double>cR);

    /**
     * @brief 计算衰减
     * @param lat1, lon1, h1 雷达的纬度、经度、高度
     * @param lat2, lon2, h2 目标的纬度、经度、高度
     * @param lata, lona, latb, lonb 左下右上两个对角点的纬度和经度
     * @param ha，hb 云底，云顶
     * @return 衰减值（dB）
     */
    double calculatePathLength(double lat1, double lon1, double h1,
                               double lat2, double lon2, double h2,
                               double lata, double lona, double latb, double lonb,
                               double ha, double hb);



private:

        //根据经度，纬度查找文件，返回对应Lred Annual Maps文件中的数据
        double findLred(const std::string& longitudeFile, const std::string& latitudeFile, const std::string& LredFile, double inputLongitude, double inputLatitude);
        // 解析文件名中的概率数字
            double parseProbability(const std::string& filename) ;
        // 计算与目标概率的差异
            double calculateDifference(double prob1, double prob2);

        //根据概率，查找目录下所有Lred文件，选择相等或者最接近的两个文件
        std::vector<std::string> getFileName(const std::string& directoryPath, double targetProbability);

        Location blhToXyz(Position p);
        //计算仰角
        double computeAngle(Location satelliteLoc, Location stationLoc, Position station);
        //判断参数合理性
        bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                Receive_polar receive_polar,double dia,double antenna_effeciency, Cloud_type cloud_type,double L,double water_tmp);

        //判断参数合理性
        bool parameter_rational1(double f,double L,double water_tmp,double lat1, double lon1,
                                 double lat2, double lon2,
                                 double lata, double lona, double latb, double lonb);
};


#endif // CLOUDATT_H
