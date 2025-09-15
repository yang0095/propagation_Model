#ifndef CISVALUEBETWEEN_H
#define CISVALUEBETWEEN_H
#include <vector>
#include <algorithm>
#include <utility> // for std::swap

class cisValueBetween
{
public:
    cisValueBetween();
    /**
     * @brief 从数组中筛选出位于指定区间内的元素
     * @param x 输入数组
     * @param a 边界点1
     * @param b 边界点2
     * @return 包含所有在区间 [min(a,b), max(a,b)] 内的元素的数组
     */
    std::vector<double> isValueBetween(const std::vector<double>& x, double a, double b) ;

    /**
     * @brief 模板版本，支持各种数值类型
     * @tparam T 数值类型
     * @param x 输入数组
     * @param a 边界点1
     * @param b 边界点2
     * @return 包含所有在区间 [min(a,b), max(a,b)] 内的元素的数组
     */
    template <typename T> std::vector<T> isValueBetween(const std::vector<T>& x, T a, T b) ;
};

#endif // CISVALUEBETWEEN_H
