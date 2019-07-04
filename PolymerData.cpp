#include "PolymerData.h"

PolymerData::PolymerData(){};

PolymerData::~PolymerData(){};

void PolymerData::initialize(double rho, long double pbh, int TunnDist, float p)
{
    Density = rho;
    PolymerBarrierHeight = pbh;
    TunnelingDistance = TunnDist;
    nu = p;
}

double PolymerData::GetDensity()
{
    return Density;
}

float PolymerData::GetPoissonsRatio()
{
    return nu;
}
