#include "cdegtorad.h"

CdegTorad::CdegTorad()
{

}

double CdegTorad::deg2rad(double deg)
{
    return deg * M_PI / 180.0;
}

double  CdegTorad::rad2deg(double rad) {
    return rad * 180.0 / M_PI;
}
