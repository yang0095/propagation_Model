#ifndef DECAYSUM_H
#define DECAYSUM_H
#include"cloudatt.h"
#include"sciatt.h"
#include"rainatt.h"
#include"fogatt.h"
#include"freespaceatt.h"
#include"hazeatt.h"
#include"snowatt.h"
#include"gasatt.h"
#include<iostream>
#include <fstream>
#include<string>

using namespace std;




class  DecaySum
{
public:
    /**
     * @brief decay_sum
     * @param fileName 计算衰减的.txt文件
     * @return  总衰减
     */
   // double decay_sum(string&fileName);

    struct StaticData{
        double rader_m;//雷达反射截面积
        double temp;//大气温度
        double humidity;//大气湿度
        double pressure;//压强
        int Rain_area_count;//雨区数量
        vector<double>Rain_rate;//每个雨区降雨率
        vector<vector<double>>Rian_area;//每个雨区左下与右上角经纬高
        int Snow_area_count;//雪区数量
        vector<int>Snow_rate;//每个雪区降雪率
        vector<vector<double>>Snow_area;//每个雪区左下与右上角经纬高
        int Fog_area_count;//雾区数量
        vector<double>Fog_rate;//每个雾区降雨率
        vector<vector<double>>Fog_area;//每个雾区左下与右上角经纬高
        int Cloud_area_count;//云区数量
        vector<double>Cloud_rate;//每个云区含水量
        vector<vector<double>>Cloud_area;//每个云区左下与右上角经纬高
        int Haze_area_count;//霾区数量
        vector<int>Haze_rate;//每个霾区污染程度
        vector<vector<double>>Haze_area;//每个霾区左下与右上角经纬高
    };
    struct DynamicData{
        double Time_Start;//仿真帧起始时间
        double Time_End;//仿真帧结束时间
        bool Data_Eff;//数据有效性
        string Decay_ID;//衰减模型ID
        bool Decay_Flag;//单双路径
        string Source_ID;//信源ID
        bool Source_State;//信源状态
        string Target_ID;//信宿ID
        bool Target_State;//信宿状态
        double Freq;     //频率
        double Source_lon;//信源经度
        double Source_lat;//信源纬度
        double Source_alt;//信源高度
        int Source_Polar;//信源极化方式
        double Rader_Gain;//雷达天线增益
        double Target_lon;//信宿经度
        double Target_lat;//信宿纬度
        double Target_alt;//信宿高度
        int Target_Polar;//信宿极化方式
        double Antenna_dia;//天线直径
        double Antenna_eta;//天线效率
    };

    struct ReturnData{
        double Time_Start;//仿真帧起始时间
        double Time_End;//仿真帧结束时间
        bool Data_Eff;//数据有效性
        string Decay_ID;//衰减模型ID
        bool Decay_Flag;//单双路径
        string Source_ID;//信源ID
        bool Source_State;//信源状态
        string Target_ID;//信宿ID
        bool Target_State;//信宿状态
        double CloudAtt;//云衰减
        double RainAtt;//雨衰减
        double FogAtt;//雾衰减
        double SnowAtt;//雪衰减
        double HazeAtt;//霾衰减
        double GasAtt;//大气衰减
        double SciAtt;//闪烁衰减
        double FreespaceAtt;//自由空间衰减
        double TotalAtt;//总衰减
    };
     StaticData Static;
//    static bool isStaticInited;
//    static DecaySum::DynamicData Dynamic;

    void receiveStaticStruct(const StaticData& userStatic);
    void receiveDynamicStruct(const DynamicData& userDynamic);


    ReturnData decay_sum(const StaticData sta,DynamicData dyn);
    ReturnData decay_sum( DynamicData dyn);
};

#endif // DECAYSUM_H
