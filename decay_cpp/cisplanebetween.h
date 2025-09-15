#ifndef CISPLANEBETWEEN_H
#define CISPLANEBETWEEN_H
#include <algorithm> // for std::min, std::max

class cisPlaneBetween
{
public:
    cisPlaneBetween();
    /**
     * @brief 判断经纬度点是否在由两个点确定的矩形区域内
     * @param lat  目标点的纬度（度）
     * @param lon  目标点的经度（度）
     * @param lat1 第一个角点的纬度（度）
     * @param lon1 第一个角点的经度（度）
     * @param lat2 第二个角点的纬度（度）
     * @param lon2 第二个角点的经度（度）
     * @return true 如果点在矩形区域内，否则 false
     */
    bool isPlaneBetween(double lat, double lon, double lat1, double lon1, double lat2, double lon2);

    /**
     * @brief 增强版：包含边界选项
     * @param includeBoundary 是否包含边界点（默认 false）
     * @return true 如果点在矩形区域内
     */
    bool isPlaneBetween(double lat, double lon, double lat1, double lon1, double lat2, double lon2, bool includeBoundary);
};

#endif // CISPLANEBETWEEN_H
