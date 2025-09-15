#include"decaysum.h"
static DecaySum::StaticData Static;
static bool isStaticInited=false;
static DecaySum::DynamicData Dynamic;


void DecaySum::receiveStaticStruct(const StaticData& userStatic){
    if(!isStaticInited){
        Static=userStatic;
        isStaticInited=true;
    }
};

void DecaySum::receiveDynamicStruct(const DynamicData& userDynamic){
    Dynamic=userDynamic;
};

DecaySum::ReturnData DecaySum::decay_sum( StaticData sta,DynamicData dyn){
    receiveStaticStruct(sta);
    receiveDynamicStruct(dyn);
    cout<<Static.pressure<<endl;

    if(Dynamic.Source_State==0 || Dynamic.Target_State==0){
        cout<<"有关机情况，无法计算衰减"<<endl;
    }


    //数据预处理
    Dynamic.Freq/=1000000000;
    Send_polar s=static_cast<Send_polar>(Dynamic.Source_Polar);
    Receive_polar r=static_cast<Receive_polar>(Dynamic.Target_Polar);
    int Polar=static_cast<Cloud_type>(Dynamic.Source_Polar);
    vector<Snow_type>snow_cast;
    for(int i=0;i<Static.Snow_rate.size();i++){
        snow_cast.push_back(static_cast<Snow_type>(Static.Snow_rate[i]));
    }
    vector<Haze_type>haze_cast;
    for(int i=0;i<Static.Haze_area_count;i++){
        haze_cast.push_back(static_cast<Haze_type>(Static.Haze_rate[i]));
    }



    CloudAtt cloud;
     double cloudatt=cloud.FUN_cloud_att(Dynamic.Freq,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                         ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,
                                         Static.Cloud_area_count,Static.Cloud_area,Static.Cloud_rate);

    RainAtt rain;
    double rainatt =rain.FUN_rain_att(Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                      ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,
                                      Static.Rain_area_count,Static.Rian_area,Static.Rain_rate,
                                      Dynamic.Freq,Polar);
    SnowAtt snow;
//    vector<Snow_type> s=static_cast<Snow_type>(Static.Snow_rate)
    double snowatt=snow.FUN_snow_att(Dynamic.Freq,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                     ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,Static.Snow_area_count,
                                     Static.Snow_area,snow_cast);
    FogAtt fog;
    double fogatt=fog.FUN_fog_att(Dynamic.Freq,Static.temp+237.15,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                  ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,Static.Fog_area_count,
                                  Static.Fog_area,Static.Fog_rate);

    HazeAtt haze;
    double hazeatt=haze.FUN_haze_att(Dynamic.Freq,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt,
                                    s,Dynamic.Antenna_eta,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,
                                     Static.Haze_area_count,Static.Haze_area,haze_cast);

    GasAtt gas;

    double gasatt=gas.gas_att(Dynamic.Freq,0.0001,Dynamic.Source_lon,Dynamic.Source_lat,Dynamic.Source_alt,s
                              ,Dynamic.Target_lon,Dynamic.Target_lat,Dynamic.Target_alt,r,Dynamic.Antenna_dia
                              ,Dynamic.Antenna_eta,Static.temp,Static.humidity,Static.pressure);


    SciAtt sci;
    double sciatt=sci.sci_att(Dynamic.Freq,0.001,Dynamic.Source_lon,Dynamic.Source_lat,Dynamic.Source_alt,s
                              ,Dynamic.Target_lon,Dynamic.Target_lat,Dynamic.Target_alt,r,Dynamic.Antenna_dia
                              ,Dynamic.Antenna_eta);




    FreeSpaceAtt freespace;
    double freeatt=freespace.freespace_att(Dynamic.Freq,0.001,Dynamic.Source_lon,Dynamic.Source_lat,Dynamic.Source_alt,s
                                           ,Dynamic.Target_lon,Dynamic.Target_lat,Dynamic.Target_alt,r,Dynamic.Antenna_dia
                                           ,Dynamic.Antenna_eta,Static.rader_m,Dynamic.Rader_Gain);


        cout<<"云衰减："<<cloudatt<<endl;
        cout<<"雨衰减："<<rainatt<<endl;
        cout<<"大气衰减："<<gasatt<<endl;
        cout<<"闪烁衰减："<<sciatt<<endl;
        cout<<"自由空间衰减："<<freeatt<<endl;
        cout<<"雾衰减："<<fogatt<<endl;
        cout<<"霾衰减："<<hazeatt<<endl;
        cout<<"雪衰减："<<snowatt<<endl;
        double totalDecay;
        if(Dynamic.Decay_Flag==1){
            cloudatt*=2;
            rainatt*=2;
            gasatt*=2;
            sciatt*=2;
            fogatt*=2;
            hazeatt*=2;
            snowatt*=2;
        }


        totalDecay=cloudatt+rainatt+gasatt+sciatt+freeatt+fogatt+hazeatt+snowatt;


        std::cout<<"总衰减"<<totalDecay<<std::endl;

        ReturnData return_data;
        return_data.Time_Start=Dynamic.Time_Start;
        return_data.Time_End=Dynamic.Time_End;
        return_data.Data_Eff=Dynamic.Data_Eff;
        return_data.Decay_ID=Dynamic.Decay_ID;
        return_data.Decay_Flag=Dynamic.Decay_Flag;
        return_data.Source_ID=Dynamic.Source_ID;
        return_data.Source_State=Dynamic.Source_State;
        return_data.Target_ID=Dynamic.Target_ID;
        return_data.Target_State=Dynamic.Target_State;
        return_data.CloudAtt=cloudatt;
        return_data.RainAtt=rainatt;
        return_data.FogAtt=fogatt;
        return_data.SnowAtt=snowatt;
        return_data.HazeAtt=hazeatt;
        return_data.GasAtt=gasatt;
        return_data.SciAtt=sciatt;
        return_data.FreespaceAtt=freeatt;
        return_data.TotalAtt=totalDecay;

        return return_data;


//    return std::round(totalDecay * 10) / 10;
}


DecaySum::ReturnData DecaySum::decay_sum( DynamicData dyn){
    receiveDynamicStruct(dyn);
    cout<<Static.pressure<<endl;

    if(Dynamic.Source_State==0 || Dynamic.Target_State==0){
        cout<<"有关机情况，无法计算衰减"<<endl;
    }


    //数据预处理
    Dynamic.Freq/=1000000000;
    Send_polar s=static_cast<Send_polar>(Dynamic.Source_Polar);
    Receive_polar r=static_cast<Receive_polar>(Dynamic.Target_Polar);
    int Polar=static_cast<Cloud_type>(Dynamic.Source_Polar);
    vector<Snow_type>snow_cast;
    for(int i=0;i<Static.Snow_rate.size();i++){
        snow_cast.push_back(static_cast<Snow_type>(Static.Snow_rate[i]));
    }
    vector<Haze_type>haze_cast;
    for(int i=0;i<Static.Haze_area_count;i++){
        haze_cast.push_back(static_cast<Haze_type>(Static.Haze_rate[i]));
    }



    CloudAtt cloud;
     double cloudatt=cloud.FUN_cloud_att(Dynamic.Freq,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                         ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,
                                         Static.Cloud_area_count,Static.Cloud_area,Static.Cloud_rate);

    RainAtt rain;
    double rainatt =rain.FUN_rain_att(Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                      ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,
                                      Static.Rain_area_count,Static.Rian_area,Static.Rain_rate,
                                      Dynamic.Freq,Polar);
    SnowAtt snow;
//    vector<Snow_type> s=static_cast<Snow_type>(Static.Snow_rate)
    double snowatt=snow.FUN_snow_att(Dynamic.Freq,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                     ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,Static.Snow_area_count,
                                     Static.Snow_area,snow_cast);
    FogAtt fog;
    double fogatt=fog.FUN_fog_att(Dynamic.Freq,Static.temp+237.15,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt
                                  ,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,Static.Fog_area_count,
                                  Static.Fog_area,Static.Fog_rate);

    HazeAtt haze;
    double hazeatt=haze.FUN_haze_att(Dynamic.Freq,Dynamic.Source_lat,Dynamic.Source_lon,Dynamic.Source_alt,
                                    s,Dynamic.Antenna_eta,Dynamic.Target_lat,Dynamic.Target_lon,Dynamic.Target_alt,
                                     Static.Haze_area_count,Static.Haze_area,haze_cast);

    GasAtt gas;

    double gasatt=gas.gas_att(Dynamic.Freq,0.0001,Dynamic.Source_lon,Dynamic.Source_lat,Dynamic.Source_alt,s
                              ,Dynamic.Target_lon,Dynamic.Target_lat,Dynamic.Target_alt,r,Dynamic.Antenna_dia
                              ,Dynamic.Antenna_eta,Static.temp,Static.humidity,Static.pressure);


    SciAtt sci;
    double sciatt=sci.sci_att(Dynamic.Freq,0.001,Dynamic.Source_lon,Dynamic.Source_lat,Dynamic.Source_alt,s
                              ,Dynamic.Target_lon,Dynamic.Target_lat,Dynamic.Target_alt,r,Dynamic.Antenna_dia
                              ,Dynamic.Antenna_eta);




    FreeSpaceAtt freespace;
    double freeatt=freespace.freespace_att(Dynamic.Freq,0.001,Dynamic.Source_lon,Dynamic.Source_lat,Dynamic.Source_alt,s
                                           ,Dynamic.Target_lon,Dynamic.Target_lat,Dynamic.Target_alt,r,Dynamic.Antenna_dia
                                           ,Dynamic.Antenna_eta,Static.rader_m,Dynamic.Rader_Gain);



        cout<<"云衰减："<<cloudatt<<endl;
        cout<<"雨衰减："<<rainatt<<endl;
        cout<<"大气衰减："<<gasatt<<endl;
        cout<<"闪烁衰减："<<sciatt<<endl;
        cout<<"自由空间衰减："<<freeatt<<endl;
        cout<<"雾衰减："<<fogatt<<endl;
        cout<<"霾衰减："<<hazeatt<<endl;
        cout<<"雪衰减："<<snowatt<<endl;
        double totalDecay;
        if(Dynamic.Decay_Flag==1){
            cloudatt*=2;
            rainatt*=2;
            gasatt*=2;
            sciatt*=2;
            fogatt*=2;
            hazeatt*=2;
            snowatt*=2;
        }


        totalDecay=cloudatt+rainatt+gasatt+sciatt+freeatt+fogatt+hazeatt+snowatt;


        std::cout<<"总衰减"<<totalDecay<<std::endl;

        ReturnData return_data;
        return_data.Time_Start=Dynamic.Time_Start;
        return_data.Time_End=Dynamic.Time_End;
        return_data.Data_Eff=Dynamic.Data_Eff;
        return_data.Decay_ID=Dynamic.Decay_ID;
        return_data.Decay_Flag=Dynamic.Decay_Flag;
        return_data.Source_ID=Dynamic.Source_ID;
        return_data.Source_State=Dynamic.Source_State;
        return_data.Target_ID=Dynamic.Target_ID;
        return_data.Target_State=Dynamic.Target_State;
        return_data.CloudAtt=cloudatt;
        return_data.RainAtt=rainatt;
        return_data.FogAtt=fogatt;
        return_data.SnowAtt=snowatt;
        return_data.HazeAtt=hazeatt;
        return_data.GasAtt=gasatt;
        return_data.SciAtt=sciatt;
        return_data.FreespaceAtt=freeatt;
        return_data.TotalAtt=totalDecay;

        return return_data;


//    return std::round(totalDecay * 10) / 10;
}














//    //（zr）除自由空间以外的其他衰减
//    std::ifstream file(filePath);
//    if (!file.is_open()) {
//        std::cout<<"打开总衰减文件失败！"<<endl;
//        return 0.0;
//    }

//    double totalDecay=0.0;
//    double cloudL,cloud_water_tmp;
//    int cloud_type_value,haze_type_value;
//    double f,PP,lat,lon,alt,satlon,satlat,satalt,dia,antenna_effeciency;
//    double ele,temp,humidity,pressure;
//    int send_polar_value,receive_polar_value;
//    int meteor_condition_value;
//    double rainrate,start_distance,end_distance,rain_height;
//    int snow_type_value;
//    double snow_start_distance,snow_end_distance;
//    int fog_meteor_condtion_value;
//    double fog_start_distance,fog_end_distance,L,fog_water_tmp;
//    Send_polar send_polar;
//    Receive_polar receive_polar;
//    Meteor_condition meteor_condition;
//    Snow_type snow_type;
//    Fog_meteor_condtion fog_meteor_condtion;
//    Haze_type haze_type;

//    Cloud_type cloud_type;
//    double antenna_gain,redar_cross;


//    //Rian Area Tpye
//        int rain_area_count;
//        vector<vector< double>>rainArea;
//        vector<double >rainType;
//    //snow area tpye
//        int snow_area_count;
//        vector<vector< double>>snowArea;
//        vector<double >snowType;
//    //cloud
//        int cloud_area_count;
//        vector<vector< double>>cloudArea;
//        vector<double >cloudType;
//    //haze
//        int haze_area_count;
//        vector<vector< double>>hazeArea;
//        vector<double >hazeType;
//    //fog
//        int fog_area_count;
//        vector<vector< double>>fogArea;
//        vector<double >fogType;

//        std::string line;
//        if(getline(file,line)){
//            std::istringstream iss(line);



//        /**
//     * 公共参数
//     * @param f 频率f（GHz） 精度1
//     * @param PP 时间概率（%） 精度0.001
//     * @param satlon 卫星经度（度）精度10-6
//     * @param satlat 卫星纬度（度）精度10-6
//     * @param satalt 飞行器高度（m）精度0.1
//     * @param send_polar  信源发射天线极化方式
//     * @param lat 站点纬度（度）精度10-6
//     * @param lon 站点经度（度）精度10-6
//     * @param alt 站点高度（单位：m) 精度0.1
//     * @param receive_polar  信宿接收天线极化方式
//     * @param dia 直径（m） 精度0.1
//     * @param antenna_effeciency 天线效率 精度0.01
//     * @param antenna_gain 天线增益精度单位(dBi)0.1  ++
//     * @param redar_cross 天线效率 精度0.01  ++
//     */

//        iss>>f>>PP>>satlon>>satlat>>satalt>>send_polar_value>>lon>>lat>>alt>>receive_polar_value>>dia>>antenna_effeciency>>antenna_gain>>redar_cross;
//        send_polar=static_cast<Send_polar>(send_polar_value);
//        receive_polar=static_cast<Receive_polar>(receive_polar_value);
//    }else{
//        std::cout<<"文件格式不正确!"<<endl;
//        file.close();
//        return 0.0;
//    }

//    if(getline(file,line)){
//         std::istringstream iss(line);
//         /**
//         * 大气参数
//         * @param ele 真实仰角ele（度） 精度0.1
//         * @param temp 温度 ℃  精度0.1
//         * @param humidity 湿度 g/m3 精度10-3
//         * @param pressure 压强 Pa  精度1
//         */
//         iss>>ele>>temp>>humidity>>pressure;
//    }else{
//        std::cout<<"文件格式不正确!"<<endl;
//        file.close();
//        return 0.0;
//    }

//    if(getline(file,line)){
//        std::istringstream iss(line);
////        /**
////        * 降雨参数
////        * @param rian_area_count雨区数量
////        * @param meteor_condition 气象条件等级
////        * @param rainrate 降雨率
////        * @param start_distance 气象条件（雨）降雨起始距离m 精度0.1
////        * @param end_distance 气象条件（雨）降雨终止距离m 精度0.1
////        * @param rain_height 雨顶高度m 精度0.1
////        */
//        /**
//        * 降雨参数（新）
//        * @param rain_area_count雨区数量
//        * @param rainArea雨区参数mx6，m为雨区数量：雨区左下角经纬高、右上角经纬高
//        * @param rainrate 降雨率mx1
//        * @param start_distance 气象条件（雨）降雨起始距离m 精度0.1
//        * @param end_distance 气象条件（雨）降雨终止距离m 精度0.1
//        * @param rain_height 雨顶高度m 精度0.1
//        */

//        //读入雨区数量rain_area_count后读取每个雨区的经纬高等信息\每个雨区的降雨类型
//        iss>>rain_area_count;
//        rainArea.reserve(rain_area_count);
//        for(int i=0;i<rain_area_count;i++){
//            vector<double> row(6);
//            for(int j=0;j<6;j++){
//                iss>>row[j];
//            }
//           rainArea.push_back(row);
//        }

//        rainType.reserve(rain_area_count);
//        for(int i=0;i<rain_area_count;i++){
//            double type;
//            iss>>type;
//            rainType.push_back(type);
//        }


////        iss>>meteor_condition_value>>rainrate>>start_distance>>end_distance>>rain_height;
////        meteor_condition=static_cast<Meteor_condition>(meteor_condition_value);
//    }else{
//        std::cout<<"文件格式不正确!"<<endl;
//        file.close();
//        return 0.0;
//    }

//    if(getline(file,line)){
//        std::istringstream iss(line);
//        /**
//        * 降雪参数
//        * @param snow_type 气象条件（雪）降雪量
//        * @param snow_start_distance   气象条件（雪）起始距离
//        * @param snow_end_distance 气象条件（雪）终止距离
//        */
//        iss>>snow_area_count;
//        snowArea.reserve(snow_area_count);
//        for(int i=0;i<snow_area_count;i++){
//            vector<double> row(6);
//            for(int j=0;j<6;j++){
//                iss>>row[j];
//            }
//           snowArea.push_back(row);
//        }

//        snowType.reserve(snow_area_count);
//        for(int i=0;i<snow_area_count;i++){
//            double type;
//            iss>>type;
//            snowType.push_back(type);
//        }
////        iss>>snow_type_value>>snow_start_distance>>snow_end_distance;
////         snow_type=static_cast<Snow_type>(snow_type_value);
//    }else{
//        std::cout<<"文件格式不正确!"<<endl;
//        file.close();
//        return 0.0;
//    }

//    if(getline(file,line)){
//        std::istringstream iss(line);
//        /**
//        * 雾参数
//     * @param fog_meteor_condtion //气象条件（雾）等级
//     * @param fog_start_distance   气象条件（雾）起始距离
//     * @param fog_end_distance 气象条件（雾）终止距离
//     * @param L 液态水含量(mm) 精度0.10.1
//     * @param fog_water_tmp 雾中液态水温度（K）精度0.1
//     */

//        iss>>fog_area_count;
//        fogArea.reserve(fog_area_count);
//        for(int i=0;i<fog_area_count;i++){
//            vector<double> row(6);
//            for(int j=0;j<6;j++){
//                iss>>row[j];
//            }
//           fogArea.push_back(row);
//        }

//        fogType.reserve(fog_area_count);
//        for(int i=0;i<fog_area_count;i++){
//            double type;
//            iss>>type;
//            fogType.push_back(type);
//        }
////        iss>>fog_meteor_condtion_value>>fog_start_distance>>fog_end_distance>>L>>fog_water_tmp;
////        fog_meteor_condtion=static_cast<Fog_meteor_condtion>(fog_meteor_condtion_value);

//    }else{
//        std::cout<<"文件格式不正确!"<<endl;
//        file.close();
//        return 0.0;
//    }

//    if(getline(file,line)){
//        std::istringstream iss(line);
//        /**
//        * 云参数+霾参数
//     * @param cloud_type 气象条件（云）类型
//     * @param cloudL 液态水含量(mm) 精度0.10.1
//     * @param cloud_water_tmp 云中液态水温度（K）精度0.1
//     * @param haze_type 气象条件（霾）
//     */

//        iss>>cloud_area_count;
//        cloudArea.reserve(cloud_area_count);
//        for(int i=0;i<cloud_area_count;i++){
//            vector<double> row(6);
//            for(int j=0;j<6;j++){
//                iss>>row[j];
//            }
//           cloudArea.push_back(row);
//        }

//        cloudType.reserve(cloud_area_count);
//        for(int i=0;i<cloud_area_count;i++){
//            double type;
//            iss>>type;
//            cloudType.push_back(type);
//        }

////        iss>>cloud_type_value>>cloudL>>cloud_water_tmp>>haze_type_value;
//        haze_type=static_cast<Haze_type>(haze_type_value);
//        cloud_type=static_cast<Cloud_type>(fog_meteor_condtion_value);
//    }else{
//        std::cout<<"文件格式不正确!"<<endl;
//        file.close();
//        return 0.0;
//    }

//    if(getline(file,line)){
//        std::istringstream iss(line);
//        iss>>haze_area_count;
//        hazeArea.reserve(haze_area_count);
//        for(int i=0;i<haze_area_count;i++){
//            vector<double> row(6);
//            for(int j=0;j<6;j++){
//                iss>>row[j];
//            }
//           hazeArea.push_back(row);
//        }

//        hazeType.reserve(haze_area_count);
//        for(int i=0;i<haze_area_count;i++){
//            double type;
//            iss>>type;
//            hazeType.push_back(type);
//        }

//    }
