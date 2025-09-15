#ifndef CEDGERECTANGLE2_H
#define CEDGERECTANGLE2_H
#include<cisvaluebetween1.h>
#include<cgreatcirclesolver.h>
#include<cisvaluebetween.h>

// 类型别名简化代码
using Point = std::pair<double, double>;
class CedgeRectangle2
{
public:
    CedgeRectangle2();

    /**
     * @brief 计算收发连线与矩形的交点（两个站点都在矩形外部的情况）
     * @param lat1, lon1 第一个站点的纬度和经度
     * @param lat2, lon2 第二个站点的纬度和经度
     * @param lata, lona 矩形左上角的纬度和经度
     * @param latb, lonb 矩形右下角的纬度和经度
     * @return 交点列表（纬度向量和经度向量）
     */
    std::pair<std::vector<double>, std::vector<double>>
    edgeRectangle2(double lat1, double lon1, double lat2, double lon2,
                   double lata, double lona, double latb, double lonb) ;

};

#endif // CEDGERECTANGLE2_H
