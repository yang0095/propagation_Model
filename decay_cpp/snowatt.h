#ifndef SNOWATT_H
#define SNOWATT_H
#include<iostream>
#include"enum.h"
#include<cmath>
#include "rainatt.h"
#include "cloudatt.h"
using namespace std;

class SnowAtt{
public:
    /**
     * 计算雪衰减
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
     * @param snow_type 气象条件（雪）降雪量
     * @param start_distance   气象条件（雪）起始距离（单位：m）
     * @param end_distance 气象条件（雪）终止距离（单位：m）
     * @return 雪衰减
     */
    double snow_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                      Receive_polar receive_polar,double dia,double antenna_effeciency, Snow_type snow_type,double start_distance,double end_distance);

    /**
     * 计算衰减率
     * @param f 频率f（单位：Hz） 精度1
     * @param snow_type 气象条件（雪）降雪量
     *
     * @return 衰减率
     */

    double specificatten(double f,Snow_type snow_type);
    /**
     * @brief 计算雪区内有效路径长度
     * @param f 频率f（单位：Hz） 精度1
     * @param snow_type 气象条件（雪）降雪量
     * @param hR 雨高（米）
     * @param elev_d 仰角（度）
     * @param lat1, lon1, h1 站点1的纬度、经度、高度
     * @param lat2, lon2, h2 站点2的纬度、经度、高度
     * @param lata, lona, latb, lonb 矩形两个对角点的纬度和经度

     * @return 雪区内有效路径长度r（米）适用于雨 、雾、雪、霾
     */
    double FUN_snow_att(double f, double lat1, double lon1, double h1,
                                  double lat2, double lon2, double h2,
                                  int snow_area_count,vector<vector<double>>snow_area,vector<Snow_type>sR);
private:
    //判断参数合理性
    bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                     Receive_polar receive_polar,double dia,double antenna_effeciency);
    bool parameter_rational1(double f,double lat1, double lon1,
                             double lat2, double lon2,
                             double lata, double lona, double latb, double lonb);
};
#endif // SNOWATT_H
