#ifndef FOGATT_H
#define FOGATT_H
#include<iostream>
#include"cloudatt.h"
#include<math.h>
#include<enum.h>
#include<cloudatt.h>
#include"rainatt.h"
using namespace  std;


class FogAtt{
public:
    /**
     * 计算雾衰减
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
     *
     * @param fog_meteor_condtion //气象条件（雾）等级
     * @param start_distance (单位：m)  气象条件（雾）起始距离
     * @param end_distance (单位：m)   气象条件（雾）终止距离
     * @param L 液态水含量(单位：mm) 精度0.10.1
     * @param water_tmp 雾中液态水温度（单位：K）精度0.1
     * @return 雾衰减
     */
    double fog_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                   Receive_polar receive_polar,double dia,double antenna_effeciency, Fog_meteor_condtion meteor_condtion,double start_distance,double end_distance,double L,double water_tmp);

    /**
     * 计算云衰减率
     * @param f 频率f（单位：Hz） 精度1
     * @param L 液态水含量(单位：mm) 精度0.10.1
     * @param water_tmp 云中液态水温度（单位：K）精度0.1
     * @param hR 雨高（米）
     * @param elev_d 仰角（度）
     * @param lat1, lon1, h1 站点1的纬度、经度、高度
     * @param lat2, lon2, h2 站点2的纬度、经度、高度
     * @param lata, lona, latb, lonb 矩形两个对角点的纬度和经度
     * @return 云衰减率
     */

    double FUN_fog_att(double f,double temp,
                       double lat1, double lon1, double h1,
                       double lat2, double lon2, double h2,
                       int fog_area_count,vector<vector<double>>fog_area, vector<double > fR);
private:
    bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                    Receive_polar receive_polar,double dia,double antenna_effeciency, Fog_meteor_condtion meteor_condtion,double start_distance,double end_distance,double L,double water_tmp);
    bool parameter_rational1(double f,double L,double water_tmp,double lat1, double lon1,
                             double lat2, double lon2,
                             double lata, double lona, double latb, double lonb);

};
#endif // FOGATT_H
