#include "ccalculateelevation.h"

CcalculateElevation::CcalculateElevation()
{

}

double CcalculateElevation::calculateElevation(double latA, double lonA, double hA, double latB, double lonB, double hB)
{
    // 将经纬度转换为弧度
    CdegTorad mydegrad;
      double latA_rad = mydegrad.deg2rad(latA);
      double lonA_rad = mydegrad.deg2rad(lonA);
      double latB_rad = mydegrad.deg2rad(latB);
      double lonB_rad = mydegrad.deg2rad(lonB);

      // 计算经度差
      double dLon = lonB_rad - lonA_rad;

      // 计算地心角γ的余弦（球面余弦公式）
      double cos_gamma = sin(latA_rad) * sin(latB_rad)
                       + cos(latA_rad) * cos(latB_rad) * cos(dLon);

      // 防止数值误差导致的越界
      if (cos_gamma > 1.0) cos_gamma = 1.0;
      if (cos_gamma < -1.0) cos_gamma = -1.0;

      // 计算两点到地心的距离
      double rA = EARTH_RADIUS + hA;
      double rB = EARTH_RADIUS + hB;

      // 计算两点间的直线距离（余弦定理）
      double d_squared = rA*rA + rB*rB - 2*rA*rB*cos_gamma;

      // 处理重合点的情况
      if (d_squared < std::numeric_limits<double>::epsilon()) {
          return std::numeric_limits<double>::quiet_NaN();
      }

      double d = std::sqrt(d_squared);

      // 计算仰角的正弦值（几何关系）
      double sin_theta = (rB * cos_gamma - rA) / d;

      // 确保正弦值在合法范围内
      if (sin_theta > 1.0) sin_theta = 1.0;
      if (sin_theta < -1.0) sin_theta = -1.0;

      // 计算仰角（弧度）
      double theta_rad = std::asin(sin_theta);

      // 转换为度数
      return mydegrad.rad2deg(theta_rad);
}
