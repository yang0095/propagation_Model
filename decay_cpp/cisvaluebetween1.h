#ifndef CISVALUEBETWEEN1_H
#define CISVALUEBETWEEN1_H
#include <utility> // for std::swap

class cisValueBetween1
{
public:
    cisValueBetween1();

    /**
     * @brief 判断数值是否在给定区间内
     * @param x 要判断的值
     * @param a 第一个边界值
     * @param b 第二个边界值
     * @return true 如果 x 在 [min(a,b), max(a,b)] 区间内，否则 false
     */
    bool isValueBetween1(double x, double a, double b) ;

    // 重载版本：使用模板支持任何可比较类型
    template <typename T>
    bool isValueBetween1_template(T x, T a, T b) ;
};

#endif // CISVALUEBETWEEN1_H
