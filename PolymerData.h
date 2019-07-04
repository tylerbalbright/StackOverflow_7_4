/* The polymer data object holds information regarding the polymer characteristics. */

#pragma once

class PolymerData
{
private:
    int TunnelingDistance;
    long double PolymerBarrierHeight;
    double Density;
    float nu;

public:
    //Constructor and Destructor
    PolymerData();
    ~PolymerData();
    
    void initialize(double rho, long double pbh, int Td, float p);
    double GetDensity();
    float GetPoissonsRatio();
};
