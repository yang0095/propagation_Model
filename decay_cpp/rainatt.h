#ifndef RAINATT_H
#define RAINATT_H
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
#include<limits>
#include"reconfig.h"
#include<enum.h>

#include"cdegtorad.h"
#include<cisplanebetween.h>
#include<cedgerectangle1.h>
#include<cedgerectangle2.h>
#include"ccalculateheight.h"
using namespace std;




class RainAtt
{
public:


    /**
     * @brief CalAttRain
     * @param f 频率f（单位：Hz） 精度1
     * @param percentage 时间概率（单位：%） 精度0.001
     * @param satlon 卫星经度（单位：度）精度10-6
     * @param satlat 卫星纬度（单位：度）精度10-6
     * @param satalt 飞行器高度（单位：m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param site_lon 站点经度（单位：度）精度10-6
     * @param site_lat 站点纬度（单位：度）精度10-6
     * @param altitude 站点海拔高度h_s（单位：）精度0.1
     * @param receive_polar  信宿接收天线极化方式
     * @param dia 直径（单位：m） 精度0.1
     * @param antenna_effeciency 天线效率 精度0.01
     * @param meteor_condition 气象条件等级
     * @param rainrate 降雨率
     * @param start_distance 气象条件（雨）降雨起始距离(单位：m)精度0.1
     * @param end_distance 气象条件（雨）降雨终止距离(单位：m) 精度0.1
     * @param rain_height 雨顶高度(单位：m) 精度0.1
     */
    double rain_att(double f,double percentage,double satlon, double satlat,double satalt, Send_polar send_polar,double site_lon, double site_lat,double altitude,
                     Receive_polar receive_polar,double dia,double antenna_effeciency, Meteor_condition meteor_condition,double rainrate,
                    double start_distance,double end_distance,double rain_height);



private:


    Location blhToXyz(Position p);

    double computeAngle(Location satelliteLoc, Location stationLoc, Position station);

    double xpd(double f,double p,double th,double tao,double ap);

    double nearestExtrapolate(const std::vector<double>&pn,const std::vector<double>&sign,double p);

    //根据经度，维度查找文件，返回对应降雨率文件中的数据
    double findValue(const std::string& longitudeFile, const std::string& latitudeFile, const std::string& nWetFile, double inputLongitude, double inputLatitude) ;

    // 实现雨衰减率的计算
    /**
     * @param theta 仰角
     * @param tao   极化方式，圆极化（45），水平极化（0），垂直极化（90）
     * @param f     频率（GHz）
     * @param r     降雨率
     * @return [0]kk [1]aa [2]att雨衰减率 dB/km
     */
    std::vector<double> specificatten(double theta, double tao, double f, double r) ;
    // 计算斜路径长度
    double GetSLength(double theta, double altitude, double RainHeight);

    //判断参数合理性
    bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                            Receive_polar receive_polar,double dia,double antenna_effeciency,double rainrate,double rain_height);

public:
    /**
     * @brief 计算雨区内有效路径长度
     * @param lat1, lon1, h1 站点1的纬度、经度、高度
     * @param lat2, lon2, h2 站点2的纬度、经度、高度
     * @param lata, lona, latb, lonb 矩形两个对角点的纬度和经度
     * @param hR 雨高（米）
     * @param elev_d 仰角（度）
     * @return 雨区内有效路径长度r（米）适用于雨 、雾、雪、霾
     */
    double calculateRainPathLength(double lat1, double lon1, double h1,
                                  double lat2, double lon2, double h2,
                                  double lata, double lona, double latb, double lonb,
                                  double hR, double elev_d);

    /**
     * @brief 计算雨衰减
     * @param lat1, lon1, h1 站点1的纬度、经度、高度
     * @param lat2, lon2, h2 站点2的纬度、经度、高度
     * @param lata, lona, latb, lonb 矩形两个对角点的纬度和经度
     * @param hR 雨高（米）
     * @param elev_d 仰角（度）
     * @param rR 降雨率（mm/h）
     * @param frequency 频率（GHz）
     * @param polarization 极化方式（0:水平，1:垂直）
     * @return 雨衰减值（dB）
     */
      double FUN_rain_att(double lat1, double lon1, double h1, double lat2, double lon2, double h2,
                          int rain_area_count,vector<vector<double>>rain_area, vector<double > rR,
                          double frequency, int polarization);
      //zy计算极化衰减量
      double GetPolar(double polar);


};

#endif // RAINATT_H
