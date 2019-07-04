/* Sub RVE data storage object. */
#include <vector>
#pragma once

class RVEdata;
class SubRVEdata
{
private:
    int ID;
    long int NumDispFail;
    std::vector<int> NumParticles;
    std::vector<std::vector<long int> > BoundaryParticles;
    int UpperX, LowerX, UpperY, LowerY, UpperZ, LowerZ;
    
public:
    //Constructor and destructor
    SubRVEdata();
    ~SubRVEdata();
    
    //Initialize instance of Sub RVE
    void initialize(int i, int Ux, int Lx, int Uy, int Ly, int Uz, int Lz, std::vector<int> NumParticles);
    
    //Place CB particles in the sub domain
    void PlaceCBParticles(int index, int NumCBFillers, RVEdata *dat, int Td, double d, float DispQuality, int Xbound);
    
    //Check RVE for corner and edge placed particles
    void CheckCornersAndEdges();
    
    //Place CNT particles in the sub domain
    void PlaceCNTParticles();
    
    //Return the boundary particle vector for the sub domain
    std::vector<std::vector<long int> > &GetBoundaryParticles();
    
    long int GetNumDispFail();
};
