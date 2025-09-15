#ifndef HAZEATT_H
#define HAZEATT_H
#include<enum.h>
#include<iostream>
#include <vector>
#include "rainatt.h"
#include "ccalculateelevation.h"
using namespace std;

class HazeAtt{
public:
    /**
     * 计算霾衰减
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
     * @param haze_type 气象条件（霾）
     * @return 霾衰减
     */
    double haze_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                      Receive_polar receive_polar,double dia,double antenna_effeciency, Haze_type haze_type);

    /**
     * @brief FUN_m_haze_att
     * @param f 频率f（单位：Hz） 精度1
     * @param satlon 卫星经度（单位：度）精度10-6
     * @param satlat 卫星纬度（单位：度）精度10-6
     * @param satalt 飞行器高度（单位：m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param ag    天线增益 dB
     * @param apd    天线物理直径 m
     * @param antenna_effeciency 天线效率 精度0.01
     * @param lon 站点经度（单位：度）精度10-6
     * @param lat 站点纬度（单位：度）精度10-6
     * @param alt 站点海拔高度h_s（单位：m）精度0.1
     * @param Rcs  目标雷达散射截面
     * @param haze_type 气象条件（霾）
     * @param LatlonArea 区域范围
     * @return 霾衰减
     */
    std::vector<double> FUN_m_haze_att(double f,double satlon, double satlat,double satalt,
                        Send_polar send_polar,double ag, double apd, double antenna_effeciency,
                        double lon, double lat,double alt,double Rcs,
                      Haze_type haze_type,std::vector< std::vector<double>> & LatlonArea);
    /**
     * @brief FUN_haze_att
     * @param f 频率f（单位：Hz） 精度1
     * @param satlon 卫星经度（单位：度）精度10-6
     * @param satlat 卫星纬度（单位：度）精度10-6
     * @param satalt 飞行器高度（单位：m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param ag    天线增益 dB
     * @param apd    天线物理直径 m
     * @param antenna_effeciency 天线效率 精度0.01
     * @param lon 站点经度（单位：度）精度10-6
     * @param lat 站点纬度（单位：度）精度10-6
     * @param alt 站点海拔高度h_s（单位：m）精度0.1
     * @param Rcs  目标雷达散射截面
     * @param haze_type 气象条件（霾）
     * @param lata, lona, latb, lonb 矩形两个对角点的纬度和经度
     * @param hR 雨高（米）
     * @return 霾衰减
     */
    double FUN_haze_att(double f,double satlat,double satlon, double satalt,
                        Send_polar send_polar, double antenna_effeciency,
                        double lat, double lon, double alt,int haze_area_count,
                        vector<vector<double>>haze_area, vector<Haze_type>hR);


private:
    //判断输入参数合理性
    bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                            Receive_polar receive_polar,double dia,double antenna_effeciency);

    bool parameter_rational1(double f,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                            double antenna_effeciency);
    double specificatten(double f,Haze_type haze_type);
};
#endif // HAZEATT_H
