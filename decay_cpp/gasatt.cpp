#include "gasatt.h"

    /**
     * @brief gas_att
     * @param f 频率f（Hz->GHz） 精度1
     * @param PP 时间概率（%） 精度0.001
     * @param satlon 卫星经度（度）精度10-6
     * @param satlat 卫星纬度（度）精度10-6
     * @param h_areo 飞行器高度（m）精度0.1
     * @param send_polar  信源发射天线极化方式
     * @param lon 站点经度（度）精度10-6
     * @param lat 站点纬度（度）精度10-6
     * @param h_s 站点海拔高度h_s（m）精度0.1
     * @param receive_polar  信宿接收天线极化方式
     * @param dia 直径（m） 精度0.1
     * @param antenna_effeciency 天线效率 精度0.01
     * @param ele 真实仰角ele（度） 精度0.1
     * @param temp 温度 ℃  精度0.1
     * @param humidity 湿度 g/m3 精度10-3
     * @param pressure 压强 Pa  精度1
     * @return 大气衰减
     */
    double GasAtt::gas_att(double f,double PP,double satlon, double satlat,double h_areo,Send_polar send_polar,double lon, double lat,double h_s,Receive_polar receive_polar,
                    double dia,double antenna_effeciency,double temp,double humidity,double pressure)
{

        CcalculateElevation mycalele;
        double ele= mycalele.calculateElevation(satlat,satlon,h_areo,lat,lon,h_s);
        if(ele <0) ele=-ele;
    //判断参数合理性
        if(parameter_rational( f, PP, satlon,  satlat, h_areo, send_polar, lon,  lat, h_s, receive_polar,dia, antenna_effeciency, ele, temp, humidity, pressure)==false){
            exit(0);
        }


    h_areo = std::min(h_areo, 20.0 * 1000.0); //只计算到20km高度以下,hlf添加
    h_s /= 1000.0; //m->km
    h_areo /= 1000.0;//m->km

    std::vector<double> h(171);
    std::vector<double> T(171);
    std::vector<double> P(171);
    std::vector<double> rou0(171);

    //计算过程
    double T0,P0,Rou0;

    T0 = temp+6.5*h_s+273.5;//c->k

    double dtemp = -6.5;
    P0 = pressure*pow(T0/(T0+dtemp*(h_s)), -34.163/dtemp);

    double H0=2.0;
    Rou0 = humidity*exp(h_s/H0);

    double gama = 0.0;
    double r0 = 6371.0;
    double e;  //定义水汽压
    int jj = 0; //分层数目
    std::vector<double> delta1(1000); //用于确定分层数目

    ele *= M_PI / 180.0; //仰角转化为弧度

    // Determine the number of layers based on aircraft altitude
    double sum_delta = h_s; // total height, cannot exceed aircraft altitude
    for (int i = 0; i < 1000; ++i) {
        delta1[i] = 0.0001 * exp(i * 1E-2);
        for (int j = 0; j < i + 1; ++j) {
            sum_delta += delta1[j];
        }
        if (sum_delta >= h_areo) {
            jj = i;
            break;
        } else {
            sum_delta = h_s;
        }
    }

    // Intermediate variables
    std::vector<double> a(1000);
    std::vector<double> alfa(1000);
    std::vector<double> beta(1000);
    std::vector<double> r(1000);
    std::vector<double> delta(1000);
    std::vector<double> n(1000);
    double Th;
    double rou_h;
    double Ph;
    double N;
    double r_abs;
    double TotalAtt = 0;

    //第一层初值
    r[0] = r0 + h_s;//第一层地球中心到第一层的高度
    beta[0] = M_PI / 2 - ele;

    for (int i = 0; i < jj; i++) {
        delta[i] = 0.0001 * exp(i * 1E-2);

        if (i == 0)//分层第一层
        {
            a[i] = -r[i] * cos(beta[i]) + 0.5 *pow(4 * pow(r[i], 2.0) * pow(cos(beta[i]), 2) + 8.0 * r[i] * delta[i] + 4.0 * pow(delta[i], 2), 0.5);
            alfa[i] = M_PI - acos(-(pow(a[i], 2) + 2 * r[i] * delta[i] + pow(delta[i], 2)) / (2 * a[i] * r[i] + 2 * a[i] * delta[i]));
            r_abs = r[i] - r0;//每层海拔高度
            Ph = MeanGlobalPres(r_abs,T0,P0);
            Th = MeanGlobalTemp(r_abs,T0);
            rou_h = MeanGlobalRou(r_abs,Rou0);
//            Th = T[i] - 6.5 * r_abs;//离地高度为r_abs时的温度,压强和水汽密度
//            Ph = P[i] * pow(T[i] / (T[i] - 6.5 * r_abs), -34.163 / 6.5);
//            rou_h = rou0[i] * exp(-r_abs / 2.0);

            e = rou_h * Th / 216.7;
            N = 77.6 / Th * (Ph + 4810.0 * e / Th);
            n[i] = 1 + N * 1E-6;
            gama = spec_gasatt(f, Th, rou_h, Ph);
            TotalAtt = a[i] * gama + TotalAtt;

        } else {
            r[i] = r[i - 1] + delta[i - 1];
            r_abs = r[i] - r0;
            Ph = MeanGlobalPres(r_abs,T0,P0);
            Th = MeanGlobalTemp(r_abs,T0);
            rou_h = MeanGlobalRou(r_abs,Rou0);
            e = rou_h * Th / 216.7;
            N = 77.6 / Th * (Ph + 4810.0 * e / Th);
            n[i] = 1 + N * 1E-6;
            beta[i] = asin(n[i - 1] / n[i] * sin(alfa[i - 1]));
            a[i] = -r[i] * cos(beta[i]) + 0.5 * pow(4.0 * pow(r[i], 2) * pow(cos(beta[i]), 2) + 8.0 * r[i] * delta[i] + 4.0 * pow(delta[i], 2), 0.5);
            alfa[i] = M_PI - acos(-(pow(a[i], 2) + 2.0 * r[i] * delta[i] + pow(delta[i], 2)) / (2.0 * a[i] * r[i] + 2.0 * a[i] * delta[i]));
            gama = spec_gasatt(f, Th, rou_h, Ph);
            if (!isnan(gama)) {
                TotalAtt = a[i] * gama + TotalAtt;
            }
        }
    }
    return std::round(TotalAtt * 10) / 10;
    }

    double GasAtt::FUN_gas_att(double f, double satlon, double satlat, double h_areo, Send_polar send_polar, double ag, double apd, double eita, double lon, double lat, double h_s, double Rcs, double temp, double humidity, double pressure)
    {
        //判断参数合理性
            if(parameter_rational1( f,  satlon,  satlat, h_areo, send_polar, lon,  lat, h_s, temp, humidity, pressure)==false){
                exit(0);
            }

        CcalculateElevation mycalele;
        f=f/pow(10,9);
        pressure/=100;//pa->hpa
        h_areo = std::min(h_areo, 20.0 * 1000.0); //只计算到20km高度以下,hlf添加
        h_s /= 1000.0; //m->km
        h_areo /= 1000.0;//m->km

        std::vector<double> h(171);
        std::vector<double> T(171);
        std::vector<double> P(171);
        std::vector<double> rou0(171);

        //计算过程
        double T0,P0,Rou0;

        T0 = temp+6.5*h_s+273.5;//c->k

        double dtemp = -6.5;
        P0 = pressure*pow(T0/(T0+dtemp*(h_s)), -34.163/dtemp);

        double H0=2.0;
        Rou0 = humidity*exp(h_s/H0);

        double gama = 0.0;
        double r0 = 6371.0;
        double e;  //定义水汽压
        int jj = 0; //分层数目
        std::vector<double> delta1(1000); //用于确定分层数目
        double ele= mycalele.calculateElevation(satlat, satlon,h_areo,lat,lon,h_s);

        ele *= M_PI / 180.0; //仰角转化为弧度

        // Determine the number of layers based on aircraft altitude
        double sum_delta = h_s; // total height, cannot exceed aircraft altitude
        for (int i = 0; i < 1000; ++i) {
            delta1[i] = 0.0001 * exp(i * 1E-2);
            for (int j = 0; j < i + 1; ++j) {
                sum_delta += delta1[j];
            }
            if (sum_delta >= h_areo) {
                jj = i;
                break;
            } else {
                sum_delta = h_s;
            }
        }

        // Intermediate variables
        std::vector<double> a(1000);
        std::vector<double> alfa(1000);
        std::vector<double> beta(1000);
        std::vector<double> r(1000);
        std::vector<double> delta(1000);
        std::vector<double> n(1000);
        double Th;
        double rou_h;
        double Ph;
        double N;
        double r_abs;
        double TotalAtt = 0;

        //第一层初值
        r[0] = r0 + h_s;//第一层地球中心到第一层的高度
        beta[0] = M_PI / 2 - ele;

        for (int i = 0; i < jj; i++) {
            delta[i] = 0.0001 * exp(i * 1E-2);

            if (i == 0)//分层第一层
            {
                a[i] = -r[i] * cos(beta[i]) + 0.5 *pow(4 * pow(r[i], 2.0) * pow(cos(beta[i]), 2) + 8.0 * r[i] * delta[i] + 4.0 * pow(delta[i], 2), 0.5);
                alfa[i] = M_PI - acos(-(pow(a[i], 2) + 2 * r[i] * delta[i] + pow(delta[i], 2)) / (2 * a[i] * r[i] + 2 * a[i] * delta[i]));
                r_abs = r[i] - r0;//每层海拔高度
                Ph = MeanGlobalPres(r_abs,T0,P0);
                Th = MeanGlobalTemp(r_abs,T0);
                rou_h = MeanGlobalRou(r_abs,Rou0);
    //            Th = T[i] - 6.5 * r_abs;//离地高度为r_abs时的温度,压强和水汽密度
    //            Ph = P[i] * pow(T[i] / (T[i] - 6.5 * r_abs), -34.163 / 6.5);
    //            rou_h = rou0[i] * exp(-r_abs / 2.0);

                e = rou_h * Th / 216.7;
                N = 77.6 / Th * (Ph + 4810.0 * e / Th);
                n[i] = 1 + N * 1E-6;
                gama = spec_gasatt(f, Th, rou_h, Ph);
                TotalAtt = a[i] * gama + TotalAtt;

            } else {
                r[i] = r[i - 1] + delta[i - 1];
                r_abs = r[i] - r0;
                Ph = MeanGlobalPres(r_abs,T0,P0);
                Th = MeanGlobalTemp(r_abs,T0);
                rou_h = MeanGlobalRou(r_abs,Rou0);
                e = rou_h * Th / 216.7;
                N = 77.6 / Th * (Ph + 4810.0 * e / Th);
                n[i] = 1 + N * 1E-6;
                beta[i] = asin(n[i - 1] / n[i] * sin(alfa[i - 1]));
                a[i] = -r[i] * cos(beta[i]) + 0.5 * pow(4.0 * pow(r[i], 2) * pow(cos(beta[i]), 2) + 8.0 * r[i] * delta[i] + 4.0 * pow(delta[i], 2), 0.5);
                alfa[i] = M_PI - acos(-(pow(a[i], 2) + 2.0 * r[i] * delta[i] + pow(delta[i], 2)) / (2.0 * a[i] * r[i] + 2.0 * a[i] * delta[i]));
                gama = spec_gasatt(f, Th, rou_h, Ph);
                if (!isnan(gama)) {
                    TotalAtt = a[i] * gama + TotalAtt;
                }
            }
        }
        TotalAtt*=2;
        return std::round(TotalAtt * 10) / 10;
    }




bool GasAtt::parameter_rational(double f,double PP,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                      Receive_polar receive_polar,double dia,double antenna_effeciency,double ele,double temp,double humidity,double pressure){

    //zr转换压强单位1hpa=100pa
    pressure*=100;

    bool flag=true;
    if(f<12||f>18){
        std::cout<<"gas_att:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"gas_att:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"gas_att:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"gas_att:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"gas_att:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"gas_att:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(dia<0||dia>5){
        std::cout<<"gas_att:天线直径="<<dia<<" 不符合范围（0~5)"<<std::endl;
        exit(0);
    }
    if(antenna_effeciency<0||antenna_effeciency>1){
        std::cout<<"gas_att:天线效率="<<antenna_effeciency<<" 不符合范围（0~1)"<<std::endl;
        exit(0);
    }
    if(ele<5||ele>80){ //判断仰角合理性
        std::cout<<"gas_att:仰角="<<ele<<"不符合范围（5~80）！"<<endl;
        exit(0);
    }
    if(temp<-80||temp>50){
        std::cout<<"gas_att:温度="<<temp<<"不符合范围（-80~50）！"<<endl;
        exit(0);
    }
    if(humidity<0||humidity>20){
        std::cout<<"gas_att:湿度="<<humidity<<"不符合范围（0~20）！"<<endl;
        exit(0);
    }
    if(pressure<0||pressure>102000){
        std::cout<<"gas_att:压强="<<pressure<<"不符合范围（0~10200）！"<<endl;
        exit(0);
    }
    return flag;
}

bool GasAtt::parameter_rational1(double f,double satlon, double satlat,double satalt, Send_polar send_polar,double lon, double lat,double alt,
                                 double temp,double humidity,double pressure){

    bool flag=true;
     f=f/pow(10,9);
    if(f<12||f>18){
        std::cout<<"gas_att:频率="<<f<<"(GHZ) 超出限定范围（12~18）!"<<endl;
        flag=false;
    }

    if(satlon<-180||satlon>180){
        std::cout<<"gas_att:信源经度="<<satlon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(satlat<-90||satlat>90){
        std::cout<<"gas_att:信源纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    double h=satlat/1000;
    if(h<0||h>36000){
        std::cout<<"gas_att:信源高度="<<h<<" 超出限定范围（0~36000）km!"<<endl;
        flag=false;
    }
    if(lon<-180||lon>180){
        std::cout<<"gas_att:信宿经度="<<lon<<" 超出限定范围（-180~180）!"<<endl;
        flag=false;
    }
    if(lat<-90||lat>90){
        std::cout<<"gas_att:信宿纬度="<<satlat<<" 超出限定范围（-90~90）!"<<endl;
        flag=false;
    }
    if(temp<-80||temp>50){
        std::cout<<"gas_att:温度="<<temp<<"不符合范围（-80~50）！"<<endl;
        exit(0);
    }
    if(humidity<0||humidity>20){
        std::cout<<"gas_att:湿度="<<humidity<<"不符合范围（0~20）！"<<endl;
        exit(0);
    }
    if(pressure<0||pressure>102000){
        std::cout<<"gas_att:压强="<<pressure<<"不符合范围（0~10200）！"<<endl;
        exit(0);
    }
    return flag;
}

double GasAtt::spec_gasatt(double freq, double T, double Rou, double P){
    double gama = 0.0;
    //函数功能:计算得到给定温湿压条件下的特征衰减;
    //输入参数:频率（GHz），温度（K），水汽密度（读数据库g/m^3）,压强（hPa）
    //输出参数:特征衰减;
    double p = 0.0;               //干空气压强(总压强为水蒸气压强和干空气压强之和即P=e+p);
    double e = 0.0;               //水蒸气压强(随高度的变化会输入)
    //	double T;              //温度(k);
    double cita = 0.0;            //温度(转换);
    //	double Rou;             //水蒸气密度;
    double SO = 0.0;          //氧气吸收线强度;
    double FO = 0.0;          //氧气吸收线线性因子;
    double SW = 0.0;          //水蒸气吸收线强度;
    double FW = 0.0;          //水蒸气吸收线线性因子;
    //double TO_att=0.0;      //氧气衰减;
    cita = 300 / T;
    //水蒸汽压强的计算；
    e = Rou * T / 216.7;
    p = P - e; //干空气压强；
    //定义变量
    double dtfo = 0.0;            //氧气线宽;
    double dtfw = 0.0;            //水蒸气线宽;
    double dalto = 0.0;           //氧气修正因子;
    double daltw = 0.0;           //水气修正因子;
    double NDf = 0.0;             //大气压力造成的氮气吸收连续带;
    double NF = 0.0;              //复折射率的虚部;

    //氧气衰减谱数据
    std::vector<double> f0 = {50.474238, 50.987749, 51.503350, 52.02141, 52.542394, 53.066907, 53.595749, 54.13, 54.671159, 55.221367, 55.783802, 56.264775, 56.363389, 56.968206, 57.612484, 58.323877, 58.44659, 59.164207, 59.590983, 60.306061, 60.434776, 61.15056, 61.800154, 62.411215, 62.48626, 62.997977, 63.568518, 64.127767, 64.678903, 65.224071, 65.764772, 66.302091, 66.83683, 67.369598, 67.900867, 68.431005, 68.960311, 118.750343, 368.498350, 424.763124, 487.249370, 715.393150, 773.839675, 834.14533};
    std::vector<double> a1 = {0.94, 2.46, 6.08, 14.14, 31.02, 64.10, 124.7, 228.00, 391.80, 631.6, 953.5, 548.9, 1344.00, 1763.00, 2141.00, 2386.0, 1457.0, 2404.0, 2112.0, 2124.0, 2461.0, 2504.0, 2298.0, 1933.0, 1517.0, 1503.0, 1087.0, 733.0, 463.0, 274.0, 153.0, 80.09, 39.46, 18.32, 8.01, 3.30, 1.28, 945.0, 67.9, 638.0, 235.0, 99.60, 671.0, 180.0};
    std::vector<double> a2 = {9.694, 8.694, 7.744, 6.844, 6.004, 5.224, 4.484, 3.814, 3.194, 2.624, 2.119, 0.015, 1.66, 1.26, 0.915, 0.626, 0.084, 0.391, 0.212, 0.212, 0.391, 0.626, 0.915, 1.26, 0.083, 1.665, 2.115, 2.62, 3.195, 3.815, 4.485, 5.225, 6.005, 6.845, 7.745, 8.695, 9.695, 0.009, 0.049, 0.044, 0.049, 0.145, 0.13, 0.147};
    std::vector<double> a3 = {8.9, 9.1, 9.4, 9.7, 9.9, 10.2, 10.5, 10.7, 11, 11.3, 11.7, 17.3, 12.0, 12.4, 12.8, 13.3, 15.2, 13.9, 14.3, 14.5, 13.6, 13.1, 12.7, 12.3, 15.4, 12.0, 11.7, 11.3, 11.0, 10.7, 10.5, 10.2, 9.9, 9.7, 9.4, 9.2, 9.0, 16.3, 19.2, 19.3, 19.2, 18.1, 18.2, 18.1};
    std::vector<double> a4 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6};
    std::vector<double> a5 = {2.4, 2.2, 1.97, 1.66, 1.36, 1.31, 2.3, 3.35, 3.74, 2.58, -1.66, 3.90, -2.97, -4.16, -6.13, -2.05, 7.48, -7.22, 7.65, -7.05, 6.97, 1.04, 5.7, 3.6, -4.98, 2.39, 1.08, -3.11, -4.21, -3.75, -2.67, -1.68, -1.69, -2.0, -2.28, -2.40, -2.5, -0.36, 0, 0, 0, 0, 0, 0};
    std::vector<double>a6 = {7.9, 7.8, 7.74, 7.64, 7.51, 7.14, 5.84, 4.31, 3.05, 3.39, 7.05, -1.13, 7.53, 7.42, 6.97, 0.51, -1.46, 2.66, -0.9, 0.81, -3.24, -0.67, -7.61, -7.77, 0.97, -7.68, -7.06, -3.32, -2.98, -4.23, -5.75, -7.0, -7.35, -7.44, -7.53, -7.6, -7.65, 0.09, 0, 0, 0, 0, 0, 0};
    //水气衰减谱数据
    std::vector<double> fw = {22.23508, 67.80396, 119.99594, 183.310091, 321.225644, 325.152919, 336.222601, 380.197372, 390.134508, 437.346667, 439.150812, 443.018295, 448.001075, 470.888947, 474.689127, 488.491133, 503.569532, 504.482692, 547.67644, 552.02096, 556.936002, 620.700807, 645.866155, 658.00528, 752.033227, 841.053973, 859.962313, 899.306675, 902.616173, 906.207325, 916.171582, 923.118427, 970.315022, 987.926764, 1780.0};
    std::vector<double> b1 = {0.113, 0.0012, 0.0008, 2.42, 0.0483, 1.499, 0.0011, 11.52, 0.0046, 0.065, 0.9218, 0.1976, 10.32, 0.3297, 1.262, 0.252, 0.039, 0.013, 9.701, 14.77, 487.4, 5.012, 0.0713, 0.3022, 239.6, 0.014, 0.1472, 0.0605, 0.0426, 0.1876, 8.34, 0.0869, 8.972, 132.1, 22300.0};
    std::vector<double> b2 = {2.143, 8.735, 8.356, 0.668, 6.181, 1.54, 9.829, 1.048, 7.35, 5.05, 3.596, 5.05, 1.405, 3.599, 2.381, 2.853, 6.733, 6.733, 0.114, 0.114, 0.159, 2.20, 8.58, 7.82, 0.396, 8.18, 7.989, 7.917, 8.432, 5.111, 1.442, 10.22, 1.92, 0.258, 0.952};
    std::vector<double>b3 = {28.11, 28.58, 29.48, 30.50, 23.03, 27.83, 26.93, 28.73, 21.52, 18.45, 21.00, 18.60, 26.32, 21.52, 23.55, 26.02, 16.12, 16.12, 26.00, 26.00, 32.10, 24.38, 18.00, 32.10, 30.6, 15.9, 30.6, 29.85, 28.65, 24.08, 26.7, 29.0, 25.5, 29.85, 176.2};
    std::vector<double> b4 = {0.69, 0.69, 0.7, 0.64, 0.67, 0.68, 0.69, 0.54, 0.63, 0.6, 0.63, 0.6, 0.66, 0.66, 0.65, 0.69, 0.61, 0.61, 0.7, 0.7, 0.69, 0.71, 0.6, 0.69, 0.68, 0.33, 0.68, 0.68, 0.7, 0.7, 0.7, 0.7, 0.64, 0.68, 0.5};
    std::vector<double> b5 = {4.8, 4.93, 4.78, 5.3, 4.69, 4.85, 4.74, 5.38, 4.81, 4.23, 4.29, 4.23, 4.84, 4.57, 4.65, 5.04, 3.89, 4.01, 4.5, 4.5, 4.11, 4.68, 4.0, 4.14, 4.09, 5.76, 4.09, 4.53, 5.1, 4.7, 4.78, 5.0, 4.94, 4.55, 30.5};
    std::vector<double>b6 = {1.0, 0.82, 0.79, 0.85, 0.54, 0.74, 0.61, 0.89, 0.55, 0.48, 0.52, 0.5, 0.67, 0.65, 0.64, 0.72, 0.43, 0.45, 1.0, 1.0, 1.0, 0.68, 0.5, 1.0, 0.84, 0.45, 0.84, 0.9, 0.95, 0.53, 0.78, 0.8, 0.67, 0.9, 1.0};
    double A1, A2, B1, B2;
    double NFo = 0;             //氧气SF乘积的和;
    double NFO;

    //计算氧气衰减效应
    //频率小于指定值时，考虑所有氧气谱线
    if (freq <= 118.750343) {
        for (int i = 0; i < 44; i++) {
            dtfo = a3[i] * 1E-4 * (p * pow(cita, 0.8 - a4[i]) + 1.1 * e * cita);//线宽
            dtfo = sqrt(pow(dtfo, 2) + 2.25E-6); //修正
            dalto = (a5[i] + a6[i] * cita) * 1E-4 * P * pow(cita, 0.8);//修正因子
            A1 = (dtfo - dalto * (f0[i] - freq)) / (pow(f0[i] - freq, 2) + pow(dtfo, 2));
            A2 = (dtfo - dalto * (f0[i] + freq)) / (pow(f0[i] + freq, 2) + pow(dtfo, 2));
            FO = freq / f0[i] * (A1 + A2);
            SO = a1[i] * pow(10, -7) * p * pow(cita, 3) * exp(a2[i] * (1 - cita)); //线强度
            NFO = SO * FO;
            NFo = NFo + NFO;
        }
    } else if (freq > 118.750343)//高于指定频率，只考虑60GHz以上谱线
    {
        for (int i = 19; i < 44; i++) {
            dtfo = a3[i] * 1E-4 * (p * pow(cita, 0.8 - a4[i]) + 1.1 * e * cita);//线宽
            dtfo = sqrt(pow(dtfo, 2) + 2.25E-6); //修正
            dalto = (a5[i] + a6[i] * cita) * 1E-4 * P * pow(cita, 0.8);//修正因子
            A1 = (dtfo - dalto * (f0[i] - freq)) / (pow(f0[i] - freq, 2) + pow(dtfo, 2));
            A2 = (dtfo - dalto * (f0[i] + freq)) / (pow(f0[i] + freq, 2) + pow(dtfo, 2));
            FO = freq / f0[i] * (A1 + A2);
            SO = a1[i] * pow(10, -7) * p * pow(cita, 3) * exp(a2[i] * (1 - cita)); //线强度
            NFO = SO * FO;
            NFo = NFo + NFO;
        }
    }
    //计算水蒸气衰减效应;
    double NFw = 0;             //水蒸气SF乘积的和;
    double NFW;
    double dtfw1;       //中间变量;
    for (int j = 0; j < 35; j++) {
        dtfw1 = b3[j] * 1E-4 * (p * pow(cita, b4[j]) + b5[j] * e * pow(cita, b6[j])); //线宽
        dtfw = 0.535 * dtfw1 + sqrt(0.217 * pow(dtfw1, 2) + 2.1316E-12 * pow(fw[j], 2) / cita);    //修正后的线宽
        daltw = 0; //修正因子;
        B1 = (dtfw - daltw * (fw[j] - freq)) / (pow(fw[j] - freq, 2) + pow(dtfw, 2));
        B2 = (dtfw - daltw * (fw[j] + freq)) / (pow(fw[j] + freq, 2) + pow(dtfw, 2));

        FW = freq / fw[j] * (B1 + B2);
        SW = b1[j] * 1E-1 * e * pow(cita, 3.5) * exp(b2[j] * (1 - cita));
        NFW = SW * FW;
        NFw = NFw + NFW;
    }

    //大气压力造成的大气吸收和deby频谱的干燥吸收连续带;
    double C1, C2, d;
    d = 5.6E-4 * P * pow(cita, 0.8);
    C1 = 6.14E-5 / (d * (1 + pow(freq / d, 2)));
    C2 = 1.4E-12 * p * pow(cita, 1.5) / (1 + 1.9E-5 * pow(freq, 1.5));
    NDf = freq * p * pow(cita, 2) * (C1 + C2);

    //求解频率依赖的复折射率的虚部;
    NF = NFo + NFw + NDf;

    //大气特征衰减;
    gama = 0.1820 * freq * NF;

    return gama;
}




/*
(zr)只用这个函数：MeanGlobalAtmo
Alt: 高度，km；
Pres
temp
Rou

(zr)说明：下列参数替换为0高度处的温湿压,通篇查找一下替换；
T0：地表温度，K;
P0: 地表压强，hPa； 已修改，表中所给单位是pa，存疑？
Rou0：地表水蒸气密度， g/m3;
*/


//(zr)高度数组可以取：0:0.5:85 km  原来程序是m；重复调用MeanGlobalAtmo
//zy回复：高度单位已是km为单位

double GasAtt::MeanGlobalTemp(double Alt,double T0)
{
    double Temp;

    double Ti, T11, T20, T32, T47, T51, T71;

    //T0 = 288.15; //(zr)输入地表温度，单位K；zy
    T11 = T0-6.5*(11-0);
    T20 = T11;
    T32 = T20+1.0*(32-20);
    T47 = T32+2.8*(47-32);
    T51 = T47;
    T71 = T51-2.8*(71-51);

    double dTemp, Hi;

    if (Alt>=0 && Alt<=11)
    {
        Ti = T0;
        dTemp = -6.5;
        Hi = 0.0;
    }
    else if (Alt <= 20)
    {
        Ti = T11;
        dTemp = 0.0;
        Hi = 11.0;
    }
    else if (Alt <= 32)
    {
        Ti = T20;
        dTemp = 1.0;
        Hi = 20.0;
    }
    else if (Alt <= 47)
    {
        Ti = T32;
        dTemp = 2.8;
        Hi = 32.0;
    }
    else if (Alt <= 51)
    {
        Ti = T47;
        dTemp = 0.0;
        Hi = 47.0;
    }
    else if (Alt <= 71)
    {
        Ti = T51;
        dTemp = -2.8;
        Hi = 51.0;
    }
    else
    {
        Ti = T71;
        dTemp = -2.0;
        Hi = 71.0;
    }

    Temp = Ti+dTemp*(Alt-Hi);

    return Temp;
}

double GasAtt::MeanGlobalPres(double Alt,double T0,double P0)
{
    double Pres;

    double Ti, T11, T20, T32, T47, T51, T71;
    double Pi, P11, P20, P32, P47, P51, P71;

    //T0 = 288.15;
    T11 = T0-6.5*(11-0);
    T20 = T11;
    T32 = T20+1.0*(32-20);
    T47 = T32+2.8*(47-32);
    T51 = T47;
    T71 = T51-2.8*(71-51);

    //P0 = 1013.25;
    P11 = P0*pow(T0/(T0-6.5*(11-0)), 34.163/(-6.5));
    P20 = P11*exp(-34.163*(20-11)/T11);
    P32 = P20*pow(T20/(T20+1.0*(32-20)), 34.163/1.0);
    P47 = P32*pow(T32/(T32+2.8*(47-32)), 34.163/2.8);
    P51 = P47*exp(-34.163*(51-47)/T47);
    P71 = P51*pow(T51/(T51-2.8*(71-51)), 34.163/(-2.8));

    double dtemp, Hi;

    if (Alt>=0 && Alt<=11)
    {
        Pi = P0;
        Ti = T0;
        dtemp = -6.5;
        Hi = 0.0;
    }
    else if (Alt <= 20)
    {
        Pi = P11;
        Ti = T11;
        dtemp = 0.0;
        Hi = 11.0;
    }
    else if (Alt <= 32)
    {
        Pi = P20;
        Ti = T20;
        dtemp = 1.0;
        Hi = 20.0;
    }
    else if (Alt <= 47)
    {
        Pi = P32;
        Ti = T32;
        dtemp = 2.8;
        Hi = 32.0;
    }
    else if (Alt <= 51)
    {
        Pi = P47;
        Ti = T47;
        dtemp = 0.0;
        Hi = 47.0;
    }
    else if (Alt <= 71)
    {
        Pi = P51;
        Ti = T51;
        dtemp = -2.8;
        Hi = 51.0;
    }
    else
    {
        Pi = P71;
        Ti = T71;
        dtemp = -2.0;
        Hi = 71.0;
    }

    if (dtemp == 0.0)
    {
        Pres = Pi*exp(-34.163*(Alt-Hi)/Ti);
    }
    else
    {
        Pres = Pi*pow(Ti/(Ti+dtemp*(Alt-Hi)), 34.163/dtemp);
    }

    return Pres;
}

double GasAtt::MeanGlobalRou(double Alt,double Rou0)
{
    double Rou;

    double H0=2.0;

    Rou = Rou0*exp(-Alt/H0);

    return Rou;
}


void GasAtt::MeanGlobalAtmo(double Alt, double *Pres, double *temp, double *Rou,double T0,double P0,double Rou0)
{
    *Pres = MeanGlobalPres(Alt,T0,P0);
    *temp = MeanGlobalTemp(Alt,T0);
    *Rou = MeanGlobalRou(Alt,Rou0);
}

