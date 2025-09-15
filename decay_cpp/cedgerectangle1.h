#ifndef CEDGERECTANGLE1_H
#define CEDGERECTANGLE1_H
#include <cmath>
#include <vector>
#include <stdexcept>
#include <utility>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include<cisvaluebetween1.h>
#include<cgreatcirclesolver.h>
#include<cisvaluebetween.h>
class CedgeRectangle1
{
public:
    CedgeRectangle1();

    // 主函数：计算收发连线与矩形的交点
    std::pair<std::vector<double>, std::vector<double>>
    edgeRectangle1(double lat1, double lon1, double lat2, double lon2,
                  double lata, double lona, double latb, double lonb);
    // 辅助打印函数
    void printResult(const std::vector<double>& lats, const std::vector<double>& lons) ;
};

#endif // CEDGERECTANGLE1_H
