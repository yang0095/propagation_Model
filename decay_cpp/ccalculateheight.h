#ifndef CCALCULATEHEIGHT_H
#define CCALCULATEHEIGHT_H
#include <cmath>
#include <vector>
#include <stdexcept>
#include<cdegtorad.h>

class CcalculateHeight
{
public:
    CcalculateHeight();

    double calculateHeight(double latA, double lonA, double hA,
                           double latB, double lonB, double hB,
                           double latC, double lonC);
};

#endif // CCALCULATEHEIGHT_H
