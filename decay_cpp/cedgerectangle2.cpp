#include "cedgerectangle2.h"

CedgeRectangle2::CedgeRectangle2()
{

}

std::pair<std::vector<double>, std::vector<double> > CedgeRectangle2::edgeRectangle2(double lat1, double lon1, double lat2, double lon2, double lata, double lona, double latb, double lonb)
{
       std::vector<double> latn;
       std::vector<double> lonn;

       cisValueBetween myvalueBt;
       cisValueBetween1 myvalueBt1;
       cGreatCircleSolver mygc;
       // 1. 处理第一条纬度边（水平边）
       if (myvalueBt1.isValueBetween1(lata, std::min(lat1, lat2), std::max(lat1, lat2))) {
//           std::cout<<
           GreatCircleSolution sol = mygc.GreatCircleSolver(lat1, lon1, lat2, lon2, lata, true);
           std::vector<double> valid_lons = myvalueBt.isValueBetween(sol.lons, std::min(lona, lonb), std::max(lona, lonb));

           for (double lon_val : valid_lons) {
               latn.push_back(lata);
               lonn.push_back(lon_val);
           }
       }

       // 2. 处理第二条纬度边（水平边）
       if (myvalueBt1.isValueBetween1(latb, std::min(lat1, lat2), std::max(lat1, lat2))) {
           GreatCircleSolution sol = mygc.GreatCircleSolver(lat1, lon1, lat2, lon2, latb, true);
           std::vector<double> valid_lons = myvalueBt.isValueBetween(sol.lons, std::min(lona, lonb), std::max(lona, lonb));

           for (double lon_val : valid_lons) {
               latn.push_back(latb);
               lonn.push_back(lon_val);
           }
       }

       // 3. 处理第一条经度边（垂直边）
       if (myvalueBt1.isValueBetween1(lona, std::min(lon1, lon2), std::max(lon1, lon2))) {
           GreatCircleSolution sol =mygc.GreatCircleSolver(lat1, lon1, lat2, lon2, lona, false);
           std::vector<double> valid_lats = myvalueBt.isValueBetween(sol.lats, std::min(lata, latb), std::max(lata, latb));

           for (double lat_val : valid_lats) {
               latn.push_back(lat_val);
               lonn.push_back(lona);
           }
       }

       // 4. 处理第二条经度边（垂直边）
       if (myvalueBt1.isValueBetween1(lonb, std::min(lon1, lon2), std::max(lon1, lon2))) {
           GreatCircleSolution sol = mygc.GreatCircleSolver(lat1, lon1, lat2, lon2, lonb, false);
           std::vector<double> valid_lats = myvalueBt.isValueBetween(sol.lats, std::min(lata, latb), std::max(lata, latb));

           for (double lat_val : valid_lats) {
               latn.push_back(lat_val);
               lonn.push_back(lonb);
           }
       }

       // 返回交点集合

       return {latn, lonn};
}
