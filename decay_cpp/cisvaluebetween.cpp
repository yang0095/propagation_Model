#include "cisvaluebetween.h"

cisValueBetween::cisValueBetween()
{

}

std::vector<double> cisValueBetween::isValueBetween(const std::vector<double> &x, double a, double b)
{
    // 确保 a <= b
       if (a > b) {
           std::swap(a, b);
       }

       // 创建结果向量
       std::vector<double> result;
       result.reserve(x.size()); // 预分配空间提高效率

       // 筛选满足条件的元素
       std::copy_if(x.begin(), x.end(), std::back_inserter(result),
                    [a, b](double value) {
                        return (value >= a) && (value <= b);
                    });

       return result;
}

template<typename T>std::vector<T> cisValueBetween::isValueBetween(const std::vector<T> &x, T a, T b)
{
    // 确保 a <= b
       if (a > b) {
           std::swap(a, b);
       }

       // 创建结果向量
       std::vector<T> result;
       result.reserve(x.size());

       // 筛选满足条件的元素
       for (const T& value : x) {
           if (value >= a && value <= b) {
               result.push_back(value);
           }
       }

       return result;
}

