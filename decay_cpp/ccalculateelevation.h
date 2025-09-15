#ifndef CCALCULATEELEVATION_H
#define CCALCULATEELEVATION_H
#include <cmath>
#include <limits>
#include <cdegtorad.h>
// 地球平均半径（米）
const double EARTH_RADIUS = 6371000.0;
class CcalculateElevation
{
public:
    CcalculateElevation();

    /**
     * 计算从点A到点B的仰角（俯仰角）
     *
     * @param latA 点A的纬度（度）
     * @param lonA 点A的经度（度）
     * @param hA 点A的海拔高度（米）
     * @param latB 点B的纬度（度）
     * @param lonB 点B的经度（度）
     * @param hB 点B的海拔高度（米）
     * @return 从点A观测点B的仰角（度），如果点重合返回NaN
     */
    double calculateElevation(
        double latA, double lonA, double hA,
        double latB, double lonB, double hB) ;
};

#endif // CCALCULATEELEVATION_H
