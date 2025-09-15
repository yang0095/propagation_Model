#include "cisplanebetween.h"

cisPlaneBetween::cisPlaneBetween()
{

}

bool cisPlaneBetween::isPlaneBetween(double lat, double lon,
                    double lat1, double lon1,
                    double lat2, double lon2) {
    // 确定纬度的最小和最大值
    double min_lat = std::min(lat1, lat2);
    double max_lat = std::max(lat1, lat2);

    // 确定经度的最小和最大值
    double min_lon = std::min(lon1, lon2);
    double max_lon = std::max(lon1, lon2);

    // 检查目标点是否在区域内
    return (lat >= min_lat) && (lat <= max_lat) &&
           (lon >= min_lon) && (lon <= max_lon);
}

/**
 * @brief 增强版：包含边界选项
 * @param includeBoundary 是否包含边界点（默认 false）
 * @return true 如果点在矩形区域内
 */
bool cisPlaneBetween::isPlaneBetween(double lat, double lon,
                    double lat1, double lon1,
                    double lat2, double lon2,
                    bool includeBoundary) {
    double min_lat = std::min(lat1, lat2);
    double max_lat = std::max(lat1, lat2);
    double min_lon = std::min(lon1, lon2);
    double max_lon = std::max(lon1, lon2);

    if (includeBoundary) {
        return (lat >= min_lat) && (lat <= max_lat) &&
               (lon >= min_lon) && (lon <= max_lon);
    } else {
        return (lat > min_lat) && (lat < max_lat) &&
               (lon > min_lon) && (lon < max_lon);
    }
}
