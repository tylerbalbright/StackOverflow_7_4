/* The CB particle class holds information about the particle size, neighbors, and location within the RVE structure. */
#include <vector>

#pragma once
class RVEdata;
class CBparticle
{
private:
    int ID;     //particle number
    int xyz[3];     //XYZ coordinates
    double radius;  //defines particle size
    std::vector<int> neighbors; //holds ID's of particles neighboring this particle
    std::vector<double> NeighborDist; //Holds dist between neigboring particles
    bool Pos;   //connected to positive terminal? (t/f)
    bool Neg;   //connected to negative terminal? (t/f)
    int MatType;    //Material type integer ID
    
public:
    //Constructor and Destructor
    CBparticle();
    ~CBparticle();
    
    //Initialize id of particle
    void Initialize(int i);
    
    //Set radius to a double value r
    void SetRadius(double r);
    
    //Get particle diameter
    int GetDiameter();
    
    //Set x, y, and z location in RVE and check closeness to terminals
    void SetXYZ(int x, int y, int z, int NumX, int TunnelDist);
    
    //Add neighbor to vectors
    void AddNeighbor(int i, double d);
    
    //Return x, y, and z coordinates
    int GetXCoordinate();
    int GetYCoordinate();
    int GetZCoordinate();
    
    //Find distance to neighbor i (pre-placement)
    double FindDistanceToNeighbor(CBparticle *i, int x, int y, int z);
    
    //Using pre-set particle info (post-placement)
    double FindDistanceToNeighbor(CBparticle *i);
    
    //Using strained particle state
    double FindDistanceToNeighbor(CBparticle *i, float eps1, float eps2, float eps3);
    
    //Return positively or negatively connected status
    bool IsPositive();
    bool IsNegative();
    
    //Set Mat Type
    void SetMatType(int i, RVEdata *dat);
    
    //Return vector of neighboring particles
    std::vector<int> &GetNeighborsVector();
    
    //Return vector of neighboring particles distances
    std::vector<double> &GetNeighborsDistVector();
    
    //Return Mat Type
    int GetMatType();
    
};
