#include "ccalculateheight.h"

CcalculateHeight::CcalculateHeight()
{

}

double CcalculateHeight::calculateHeight(double latA, double lonA, double hA, double latB, double lonB, double hB, double latC, double lonC)
{
    const double R = 6371000.0;  // 地球半径(米)
    CdegTorad mydegrad;

       // 转换角度到弧度
       double latA_rad = mydegrad.deg2rad(latA);
       double lonA_rad = mydegrad.deg2rad(lonA);
       double latB_rad = mydegrad.deg2rad(latB);
       double lonB_rad = mydegrad.deg2rad(lonB);
       double latC_rad = mydegrad.deg2rad(latC);
       double lonC_rad = mydegrad.deg2rad(lonC);

       // 计算点A的直角坐标
       double xA = (R + hA) * cos(latA_rad) * cos(lonA_rad);
       double yA = (R + hA) * cos(latA_rad) * sin(lonA_rad);
       double zA = (R + hA) * sin(latA_rad);

       // 计算点B的直角坐标
       double xB = (R + hB) * cos(latB_rad) * cos(lonB_rad);
       double yB = (R + hB) * cos(latB_rad) * sin(lonB_rad);
       double zB = (R + hB) * sin(latB_rad);

       // 计算方向向量
       double dx = xB - xA;
       double dy = yB - yA;
       double dz = zB - zA;

       // 计算点C的单位方向向量分量
       double a = cos(latC_rad) * cos(lonC_rad);
       double b = cos(latC_rad) * sin(lonC_rad);
       double c = sin(latC_rad);

       // 构建线性方程组：最小二乘求解
       // 方程组形式: [a, -dx; b, -dy; c, -dz] * [k; t] = [xA; yA; zA]
       // 正规方程: (A^T A) X = A^T B

       // 计算 A^T A (2x2矩阵)
       std::vector<std::vector<double>> ATA(2, std::vector<double>(2, 0.0));
       ATA[0][0] = a*a + b*b + c*c;
       ATA[0][1] = ATA[1][0] = -(a*dx + b*dy + c*dz);
       ATA[1][1] = dx*dx + dy*dy + dz*dz;

       // 计算 A^T B (2x1向量)
       std::vector<double> ATB(2, 0.0);
       ATB[0] = a*xA + b*yA + c*zA;
       ATB[1] = -(dx*xA + dy*yA + dz*zA);

       // 计算行列式
       double det = ATA[0][0] * ATA[1][1] - ATA[0][1] * ATA[1][0];

       // 检查矩阵是否奇异
       if (std::abs(det) < 1e-12) {
           throw std::runtime_error("Matrix is singular, cannot solve the equation.");
       }

       // 计算逆矩阵: inv(ATA)
       std::vector<std::vector<double>> invATA(2, std::vector<double>(2, 0.0));
       invATA[0][0] = ATA[1][1] / det;
       invATA[1][1] = ATA[0][0] / det;
       invATA[0][1] = -ATA[0][1] / det;
       invATA[1][0] = -ATA[1][0] / det;

       // 求解 X = inv(ATA) * ATB
       double k = invATA[0][0]*ATB[0] + invATA[0][1]*ATB[1];
       double t = invATA[1][0]*ATB[0] + invATA[1][1]*ATB[1];

       // 计算高度
       double hC = k - R;

       // 计算残差
       double res1 = a*k - dx*t - xA;
       double res2 = b*k - dy*t - yA;
       double res3 = c*k - dz*t - zA;
       double residual = std::sqrt(res1*res1 + res2*res2 + res3*res3);

       if (residual > 1e-6) {
           throw std::runtime_error("Solution inaccurate, residual: " + std::to_string(residual));
       }

       return hC;
}
