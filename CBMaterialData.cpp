#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include "CBMaterialData.h"

CBMaterialData::CBMaterialData(){};

CBMaterialData::~CBMaterialData(){};

void CBMaterialData::initialize(int i, double rho, int l, int u, int ps)
{
    ID = i;
    Density = rho;
    diamLowerBound = 30;
    diamUpperBound = 30;
    PSoftness = ps;
}

double CBMaterialData::GetDensity()
{
    return Density;
}

long double CBMaterialData::GetAvgParticleWeight()
{
    long double M1p = (4/3)*M_PI*(pow((diamLowerBound+diamUpperBound)*pow(10,-9)/4,3))*Density; //Mass of 1 particle based on a spherical shape and bulk density
    return M1p;
}

void CBMaterialData::SetNumParticles(int i)
{
    NumberOfParticlesInRVE = i;
}

int CBMaterialData::GetNumParticles()
{
    return NumberOfParticlesInRVE;
}

int CBMaterialData::GetPenetrationAllowance()
{
    return PSoftness;
}

int CBMaterialData::GetRandomRadius()
{
    double d = 0;
    std::uniform_int_distribution<unsigned int> distribution(diamLowerBound, diamUpperBound);
    srand(time(NULL));
    std::mt19937_64 gen(rand());
    d = distribution(gen);
    
    return round(d/2.0);
}

int CBMaterialData::GetMaxDiameter()
{
    return diamUpperBound;
}
