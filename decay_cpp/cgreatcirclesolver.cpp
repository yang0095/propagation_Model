#include "cgreatcirclesolver.h"

// 表示大圆解的结构体
using namespace std;
cGreatCircleSolver::cGreatCircleSolver()
{

}

double cGreatCircleSolver::normalizeLongitude(double lon_deg)
{
    // 减少周期数
       lon_deg = fmod(lon_deg, 360.0);

       // 调整到 [-180, 180]
       if (lon_deg > 180.0) {
           lon_deg -= 360.0;
       } else if (lon_deg < -180.0) {
           lon_deg += 360.0;
       }

       return lon_deg;
}

GreatCircleSolution cGreatCircleSolver::GreatCircleSolver(double latA, double lonA, double latB, double lonB, double givenValue, bool isLatitudeGiven)
{
    // 容差值
      constexpr double eps = 1e-10;
      GreatCircleSolution solution;
      CdegTorad mydegtorad;

      // 转换为弧度
      double phiA = mydegtorad.deg2rad(latA);
      double lambdaA = mydegtorad.deg2rad(lonA);
      double phiB = mydegtorad.deg2rad(latB);
      double lambdaB = mydegtorad.deg2rad(lonB);

      // 转换为笛卡尔坐标（单位球）
      double xA = cos(phiA) * cos(lambdaA);
      double yA = cos(phiA) * sin(lambdaA);
      double zA = sin(phiA);

      double xB = cos(phiB) * cos(lambdaB);
      double yB = cos(phiB) * sin(lambdaB);
      double zB = sin(phiB);

      // 计算大圆平面法向量 (A × B)
      double Nx = yA * zB - zA * yB;
      double Ny = zA * xB - xA * zB;
      double Nz = xA * yB - yA * xB;

      // 检查法向量模长
      double N_norm = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
      if (N_norm < eps) {
          throw std::runtime_error("点A和点B重合或对跖，大圆未定义");
      }

      // 规范法向量
      Nx /= N_norm;
      Ny /= N_norm;
      Nz /= N_norm;

      if (isLatitudeGiven) {
          // 已知纬度求经度
          double phiP = mydegtorad.deg2rad(givenValue);

          // 处理极点情况
          if (fabs(fabs(phiP) - M_PI_2) < eps) {
              if (fabs(Nz) < eps) {
                  // 大圆经过极点
                  solution.lats.push_back(givenValue);
                  solution.lons.push_back(std::numeric_limits<double>::quiet_NaN());
                  std::cerr << "警告：大圆经过极点，经度未定义" << std::endl;
                  return solution;
              } else {
                  // 大圆不经过该极点
                  throw std::runtime_error("无解：大圆不经过该极点");
              }
          }

          // 计算方程参数
          double D = -Nz * tan(phiP);
          double R = sqrt(Nx*Nx + Ny*Ny);

          if (R < eps) {
              if (fabs(D) < eps) {
                  // 大圆是赤道
                  solution.lats.push_back(givenValue);
                  solution.lons.push_back(0.0); // 任意经度
                  std::cerr << "警告：大圆是赤道，经度任意" << std::endl;
                  return solution;
              } else {
                  throw std::runtime_error("无解：纬度不在大圆上");
              }
          }

          // 计算两个可能的经度解
          double alpha = atan2(Ny, Nx);
          double ratio = D / R;

          if (fabs(ratio) > 1.0 + eps) {
              throw std::runtime_error("无解：纬度不在大圆上");
          }

          // 限制ratio在[-1,1]范围内
          if(ratio>1.0)
          {
              ratio=1.0;
          }else if (ratio<-1.0)
          {
              ratio=-1.0;
          }
          double delta = acos(ratio);

          double lonP_rad1 = alpha + delta;
          double lonP_rad2 = alpha - delta;

          // 转换为度并规范经度
          double lon1 = normalizeLongitude(mydegtorad.rad2deg(lonP_rad1));
          double lon2 = normalizeLongitude(mydegtorad.rad2deg(lonP_rad2));

          // 存储结果
          solution.lats = {givenValue, givenValue};
          solution.lons = {lon1, lon2};
      } else {
          // 已知经度求纬度
          double lambdaP = mydegtorad.deg2rad(givenValue);

          // 计算K值
          double K = Nx * cos(lambdaP) + Ny * sin(lambdaP);

          if (fabs(Nz) < eps) {
              if (fabs(K) < eps) {
                  // 整个经线在大圆上（极点情况）
                  solution.lats = {90.0, -90.0};
                  solution.lons = {givenValue, givenValue};
                  std::cerr << "警告：整个经线在大圆上，纬度任意" << std::endl;
                  return solution;
              } else {
                  // 无解
                  throw std::runtime_error("无解：经度不在大圆上");
              }
          }

          // 计算纬度
          double tan_phi = -K / Nz;
          double phiP = atan(tan_phi);
          double latP = mydegtorad.rad2deg(phiP);

          // 存储单解
          solution.lats = {latP};
          solution.lons = {givenValue};
      }

      return solution;
}

void cGreatCircleSolver::printSolution(const GreatCircleSolution &solution)
{
    if (!solution.hasSolution()) {
           std::cout << "无解" << std::endl;
           return;
       }

       for (size_t i = 0; i < solution.numSolutions(); ++i) {
           std::cout << "解 " << (i+1) << ": ";
           std::cout << "纬度 = " << std::setw(10) << solution.lats[i];

           if (std::isnan(solution.lons[i])) {
               std::cout << ", 经度 = NaN";
           } else {
               std::cout << ", 经度 = " << std::setw(10) << solution.lons[i];
           }
           std::cout << std::endl;
       }
}
