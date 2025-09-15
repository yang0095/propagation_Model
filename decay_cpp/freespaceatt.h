#ifndef FREESPACEATT_H
#define FREESPACEATT_H
#include<iostream>
#include<math.h>
#include<enum.h>
using namespace std;

class  FreeSpaceAtt{
  public:
    /**
     * 计算自由空间传播衰减
     * @param f 频率f（Hz） 精度1
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
     * @return 自由空间传播衰减
     */
    double freespace_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                      Receive_polar receive_polar,double dia,double antenna_effeciency,double radar,double gain);

    /**
     * 计算自由空间传播衰减
     * @param f 频率f（Hz） 精度1
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
     * @param sigma 目标散射截面㎡
     * @return 自由空间传播衰减
     */

    double FUN_freespace_att(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                             Receive_polar receive_polar,double dia,double antenna_effeciency,double sigma = 1.0);


  private:
    //判断参数和理性
    bool parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                            Receive_polar receive_polar,double dia,double antenna_effeciency);
  public:
      //计算球面坐标距离，alt1单位为km
    double haversine(double lat1,double lon1,double alt1,double lat2,double lon2,double alt2);
};
#endif // FREESPACEATT_H
