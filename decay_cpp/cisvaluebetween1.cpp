#include "cisvaluebetween1.h"
using namespace std;
cisValueBetween1::cisValueBetween1()
{

}

bool cisValueBetween1::isValueBetween1(double x, double a, double b)
{
    // 确保 a <= b
       if (a > b) {
           std::swap(a, b);
       }

       // 判断 x 是否在区间内
       return (x >= a) && (x <= b);
}

template<typename T>
bool cisValueBetween1::isValueBetween1_template(T x, T a, T b)
{
    if (a > b) {
           std::swap(a, b);
       }
       return (x >= a) && (x <= b);
}
