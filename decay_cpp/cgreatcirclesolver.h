#ifndef CGREATCIRCLESOLVER_H
#define CGREATCIRCLESOLVER_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include<cdegtorad.h>
#include<utility>


struct GreatCircleSolution {
    std::vector<double> lats; // 纬度值 (度)
    std::vector<double> lons; // 经度值 (度)

    // 检查是否有解
    bool hasSolution() const {
        return !lats.empty() && !lons.empty();
    }

    // 获取解的数量
    size_t numSolutions() const {
        return std::min(lats.size(), lons.size());
    }
};
class cGreatCircleSolver
{
public:
    cGreatCircleSolver();

    // 规范经度到 [-180, 180] 范围
    double normalizeLongitude(double lon_deg);

    /**
     * @brief 大圆求解器：已知一点坐标求另一个坐标
     * @param latA 点A纬度 (度)
     * @param lonA 点A经度 (度)
     * @param latB 点B纬度 (度)
     * @param lonB 点B经度 (度)
     * @param givenValue 已知的纬度或经度 (度)
     * @param isLatitudeGiven 是否为已知纬度
     * @return GreatCircleSolution 包含解的结构体
     */
    GreatCircleSolution GreatCircleSolver(double latA, double lonA,
                                         double latB, double lonB,
                                         double givenValue,
                                         bool isLatitudeGiven);
    // 打印解决方案
    void printSolution(const GreatCircleSolution& solution);
};

#endif // CGREATCIRCLESOLVER_H
