#ifndef ENUM_H
#define ENUM_H

enum class Meteor_condition{//气象（降雨）条件等级
    Min=0,
    Medium=1,
    large=2
};

enum Cloud_type{  //气象条件（云）类型
    Cumulus=0,
    Stratus=1,
    Stratocumulus=2,
    Altostratus=3,
    Plumocumulus=4,
    Cirrus=5
};

enum Fog_meteor_condtion{ //气象条件（雾）等级
    Light=0,
    Strong=1
};

enum class Snow_type{  //气象条件（雪）类型
    None=-1,
    Min=0,
    Medium=1,
    large=2
};

enum class Send_polar{   //信源发射天线极化方式
    HH=1,
    HV=2,
    VH=3,
    VV=4
};

enum class Receive_polar{ //信宿接收天线极化方式
    HH=1,
    HV=2,
    VH=3,
    VV=4
};

enum  Haze_type{  //气象条件（霾）类型
   Slight=-1,
    Mild=0,
    Moderate=1,
    Severe=2
};

class  Position {
public:
    double lat, lon, alt;
    void setLat(double latitude) { lat = latitude; }
    void setLon(double longitude) { lon = longitude; }
    void setAlt(double altitude) { alt = altitude; }
};

class Location {
public:
    double x, y, z;
};


#endif // ENUM_H
