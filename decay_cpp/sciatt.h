#ifndef SCIATT_H
#define SCIATT_H
#include"reconfig.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include<iomanip>
#include<queue>
#include<functional>
#include<enum.h>
#include "rainatt.h"
using namespace std;


class  SciAtt
{
private:
    //根据经度，维度查找文件，返回对应n_wet文件中的数据
    double findNWet(const std::string& longitudeFile, const std::string& latitudeFile, const std::string& nWetFile, double inputLongitude, double inputLatitude);
    //格式化文件名
    std::string getFileName(double probability);


    Location blhToXyz(Position p);

    double computeAngle(Location satelliteLoc, Location stationLoc, Position station);

    //判断参数合理性
    bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                    Receive_polar receive_polar,double dia,double antenna_effeciency);

    /**
    * @param lat1, lon1, h1 站点1的纬度、经度、高度
    * @param lat2, lon2, h2 站点2的纬度、经度、高度
    * @param lata, lona, latb, lonb 矩形两个对角点的纬度和经度
    */
    bool parameter_rational1(double lat1, double lon1,
                             double lat2, double lon2,
                             double lata, double lona, double latb, double lonb);
 public:
    /**
     * 计算闪烁衰减
     * @param freq_GHz 频率f（单位：Hz） 精度1
     * @param PP 时间概率（单位：%） 精度0.001
     * @param lonSatellite 卫星经度（单位：度）精度10-6
     * @param latSatellite 卫星纬度（单位：度）精度10-6
     * @param satalt 飞行器高度（单位：m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param lonStation 站点经度（单位：度）精度10-6
     *
     *
     * @param latStation 站点纬度（单位：度）精度10-6
     * @param altStation 站点高度（单位：m) 精度0.1
     * @param receive_polar  信宿接收天线极化方式
     * @param D 直径（单位：m） 精度0.1
     * @param eita 天线效率 精度0.01
     * @return
     */
    double sci_att( double f,double PP,double lonSatellite,double latSatellite,double satalt,
                     Send_polar send_polar,double lonStation, double latStation,double altStation, Receive_polar receive_polar,double D,double eita);


    /**
     * @brief 计算闪烁衰减
     * @param lon1,lat1,  h1 信源的经度、纬度、高度
     * @param lon2,lat2,  h2 信宿的经度、纬度、高度
     * @param lata, lona, latb, lonb 矩形两个对角点的纬度和经度
     * @param temperature 温度摄氏度
     * @param elev_d 仰角（度）
     * @param rou 湿度
     * @param hR 顶高（米） 默认1000米
     * @param PP 时间概率（单位：%） 精度0.001  默认50%
     * @param Pressure 压强（hPa） 默认1013
     * @return 衰减值
     */
    double FUN_sci_att(double lat1, double lon1, double h1,
                       double lat2, double lon2, double h2,
                       double lata, double lona, double latb, double lonb,
                       double temperature, double elev_d, double rou,double hR = 1000,
                        double PP = 50,double Pressure = 1013);
};

#endif // SCIATT_H
