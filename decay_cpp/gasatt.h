#ifndef GASATT_H
#define GASATT_H
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
#include<enum.h>
#include <cstring>
#include "ccalculateelevation.h"
using namespace std;

class GasAtt
{
public:
    /**
     * @brief gas_att
     * @param f 频率f（单位：Hz） 精度1
     * @param PP 时间概率（单位：%） 精度0.001
     * @param satlon 卫星经度（单位：度）精度10-6
     * @param satlat 卫星纬度（单位：度）精度10-6
     * @param h_areo 飞行器高度（单位：m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param lon 站点经度（单位：度）精度10-6
     * @param lat 站点纬度（单位：度）精度10-6
     * @param h_s 站点海拔高度h_s（单位：m）精度0.1
     * @param receive_polar  信宿接收天线极化方式
     * @param dia 直径（单位：m） 精度0.1
     * @param antenna_effeciency 天线效率 精度0.01
     * @param ele 真实仰角ele（单位：度） 精度0.1
     * @param temp 温度 (单位：℃)  精度0.1
     * @param humidity 湿度 (单位：g/m3) 精度10-3
     * @param pressure 压强 (单位：Pa)  精度1
     * @return 大气衰减
     */
    double gas_att(double f,double PP,double satlon, double satlat,double h_areo,Send_polar send_polar,double lon, double lat,double h_s,Receive_polar receive_polar,
                    double dia,double antenna_effeciency,double temp,double humidity,double pressure);

    /**
     * @brief FUN_gas_att
     * @param f 频率f（单位：Hz） 精度1
     * @param satlon 卫星经度（单位：度）精度10-6
     * @param satlat 卫星纬度（单位：度）精度10-6
     * @param h_areo 飞行器高度（单位：m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param ag    天线增益 dB
     * @param apd    天线物理直径 m
     * @param eita 天线效率 精度0.01
     * @param lon 站点经度（单位：度）精度10-6
     * @param lat 站点纬度（单位：度）精度10-6
     * @param h_s 站点海拔高度h_s（单位：m）精度0.1
     * @param Rcs  目标雷达散射截面
     * @param temp 温度 (单位：℃)  精度0.1
     * @param humidity 湿度 (单位：g/m3) 精度10-3
     * @param pressure 压强 (单位：Pa)  精度1
     * @return 大气衰减
     */
    double FUN_gas_att(double f, double satlon, double satlat,double h_areo,Send_polar send_polar,
                       double ag, double apd, double eita,
                       double lon, double lat, double h_s,double Rcs,
                       double temp,double humidity,double pressure);



private:
    double spec_gasatt(double freq, double T, double Rou, double P);

    void MeanGlobalAtmo(double Alt, double *Pres, double *temp, double *Rou,double T0,double P0,double Rou0);

    double MeanGlobalTemp(double Alt,double T0);

    double MeanGlobalPres(double Alt,double T0,double P0);

    double MeanGlobalRou(double Alt,double Rou0);

    //判断输入参数合理性
    bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                            Receive_polar receive_polar,double dia,double antenna_effeciency,double ele,double temp,double humidity,double pressure);
    bool parameter_rational1(double f, double satlon, double satlat, double satalt, Send_polar send_polar, double lon, double lat, double alt, double temp, double humidity, double pressure);
};

#endif // GASATT_H
