#include "cedgerectangle1.h"

CedgeRectangle1::CedgeRectangle1()
{

}

std::pair<std::vector<double>, std::vector<double> > CedgeRectangle1::edgeRectangle1(double lat1, double lon1, double lat2, double lon2, double lata, double lona, double latb, double lonb)
{
    cisValueBetween1 myvalueB1;
    cGreatCircleSolver mygreatC;
    cisValueBetween myvalueBeet;
    std::vector<double> result_lats;
       std::vector<double> result_lons;

       // 确定矩形边界
       double min_lat = std::min(lata, latb);
       double max_lat = std::max(lata, latb);
       double min_lon = std::min(lona, lonb);
       double max_lon = std::max(lona, lonb);

       // 1. 尝试上下边界（固定纬度）
       if (myvalueB1.isValueBetween1(lata, lat1, lat2)) {
           GreatCircleSolution sol = mygreatC.GreatCircleSolver(lat1, lon1, lat2, lon2, lata, true);
           if (!sol.lons.empty()) {
               std::vector<double> valid_lons = myvalueBeet.isValueBetween(sol.lons, min_lon, max_lon);
               if (!valid_lons.empty()) {
                   result_lats.push_back(lata);
                   result_lons.push_back(valid_lons[0]); // 取第一个有效解
               }
           }
       }

       if (!result_lats.empty()) {
           return {result_lats, result_lons};
       }

       if (myvalueB1.isValueBetween1(latb, lat1, lat2)) {
           GreatCircleSolution sol = mygreatC.GreatCircleSolver(lat1, lon1, lat2, lon2, latb, true);
           if (!sol.lons.empty()) {
               std::vector<double> valid_lons = myvalueBeet.isValueBetween(sol.lons, min_lon, max_lon);
               if (!valid_lons.empty()) {
                   result_lats.push_back(latb);
                   result_lons.push_back(valid_lons[0]);
               }
           }
       }

       if (!result_lats.empty()) {
           return {result_lats, result_lons};
       }

       // 2. 尝试左右边界（固定经度）
       if (myvalueB1.isValueBetween1(lona, lon1, lon2)) {
           GreatCircleSolution sol = mygreatC.GreatCircleSolver(lat1, lon1, lat2, lon2, lona, false);
           if (!sol.lats.empty()) {
               std::vector<double> valid_lats =myvalueBeet.isValueBetween(sol.lats, min_lat, max_lat);
               if (!valid_lats.empty()) {
                   result_lats.push_back(valid_lats[0]);
                   result_lons.push_back(lona);
               }
           }
       }

       if (!result_lats.empty()) {
           return {result_lats, result_lons};
       }

       if (myvalueB1.isValueBetween1(lonb, lon1, lon2)) {
           GreatCircleSolution sol = mygreatC.GreatCircleSolver(lat1, lon1, lat2, lon2, lonb, false);
           if (!sol.lats.empty()) {
               std::vector<double> valid_lats = myvalueBeet.isValueBetween(sol.lats, min_lat, max_lat);
               if (!valid_lats.empty()) {
                   result_lats.push_back(valid_lats[0]);
                   result_lons.push_back(lonb);
               }
           }
       }

       return {result_lats, result_lons};
}

void CedgeRectangle1::printResult(const std::vector<double> &lats, const std::vector<double> &lons)
{
    if (lats.empty()) {
           std::cout << "无交点" << std::endl;
           return;
       }

       for (size_t i = 0; i < lats.size(); ++i) {
           std::cout << "交点 " << (i+1) << ": ";
           std::cout << "纬度 = " << std::fixed << std::setprecision(6) << lats[i];
           std::cout << ", 经度 = " << lons[i] << std::endl;
       }
}
