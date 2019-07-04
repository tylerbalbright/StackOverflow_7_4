#include <random>
#include <iostream>
#include <algorithm>
#include "RVEdata.h"
#include "SubRVEdata.h"

SubRVEdata::SubRVEdata()
{
    BoundaryParticles.resize(26);
};

SubRVEdata::~SubRVEdata(){};

void SubRVEdata::initialize(int i, int Lx, int Ux, int Ly, int Uy, int Lz, int Uz, std::vector<int> NP)
{
    for (int i = 0; i < NP.size(); i++)
    {
        NumParticles.push_back(NP[i]);
    }
    UpperX = Ux;
    LowerX = Lx;
    UpperY = Uy;
    LowerY = Ly;
    UpperZ = Uz;
    LowerZ = Lz;
    ID = i;
}

void SubRVEdata::PlaceCBParticles(int index, int NumCBFillers, RVEdata *dat, int Td, double d, float DispQuality, int Xbound)
{
    NumDispFail = 0;
    std::uniform_int_distribution<unsigned int> distributionx(LowerX, UpperX);
    std::uniform_int_distribution<unsigned int> distributiony(LowerY, UpperY);
    std::uniform_int_distribution<unsigned int> distributionz(LowerZ, UpperZ);
    
    int TotalParticles = 0;
    for (int i = 0; i < NumParticles.size(); i++)
    {
        TotalParticles += NumParticles[i];
    }
    //calculate starting index of this RVE in terms of global particles
    long int ParticleIndexOffset = (index * TotalParticles);
    
    //Find total number of CB particles to place
    TotalParticles = 0;
    for (int i = 0; i < NumCBFillers; i++)
    {
        TotalParticles += NumParticles[i];
    }
    
    srand(time(NULL));
    std::mt19937_64 gen(rand());

    //Temporary Neighbor data
    std::vector<int> neighIndex;
    std::vector<double> neighDist;

//Loop that makes sure the particles aren't inside of each other to a certain degree
    unsigned int xi=0, yi=0, zi=0;
    int SetType = 0;
    double dist;
    
    for (unsigned int i = ParticleIndexOffset; i<(ParticleIndexOffset+TotalParticles); i++)
    {
        if (i >= NumParticles[SetType]+ParticleIndexOffset)
        {
            SetType += 1;
            std::cout << "This shouldn't be happening" << std::endl;
        }
        dat->GetCBParticle(dat->GetParticle(i)->GetTypeID())->SetMatType(SetType, dat);
                
        int iPSoftness = dat->GetCBMatType(dat->GetCBParticle(dat->GetParticle(i)->GetTypeID())->GetMatType())->GetPenetrationAllowance();
        int dpi = dat->GetCBParticle(dat->GetParticle(i)->GetTypeID())->GetDiameter();
        int dpj;
        
        int t = 1;
        while (t > 0)
        {
            xi = distributionx(gen);
            yi = distributiony(gen);
            zi = distributionz(gen);
            
            int xx = xi;
            int yy = yi;
            int zz = zi;
            
            t = 0;                          //Variable to be changed if tests deny particle coords
            double smallestdist = 100;        // initializing a largest dist identifier
            
            // For loop that checks the particle coords of particles of same type  already placed
            for (unsigned int j = ParticleIndexOffset; j<i; j++)
            {
                int xj, yj, zj;
                xj = dat->GetCBParticle(dat->GetParticle(j)->GetTypeID())->GetXCoordinate();
                yj = dat->GetCBParticle(dat->GetParticle(j)->GetTypeID())->GetYCoordinate();
                zj = dat->GetCBParticle(dat->GetParticle(j)->GetTypeID())->GetZCoordinate();
                
                int jPSoftness = dat->GetCBMatType(dat->GetCBParticle(dat->GetParticle(j)->GetTypeID())->GetMatType())->GetPenetrationAllowance();
                
                //Define the most restrictive penetration allowance
                double PS;
                if (iPSoftness < jPSoftness)
                {
                    PS = iPSoftness;
                }
                else
                {
                    PS = jPSoftness;
                }
                
                dpj = dat->GetCBParticle(dat->GetParticle(j)->GetTypeID())->GetDiameter();
                
                // Define distance from particle(j) to particle(i)
                dist = dat->GetCBParticle(dat->GetParticle(i)->GetTypeID())->FindDistanceToNeighbor(dat->GetCBParticle(dat->GetParticle(j)->GetTypeID()), xx, yy, zz);
                
                // checking distance w.r.t the particle sizes and tunneling allowance
                if (dist < ((dpi+dpj)/2.0) + Td || std::abs(dist - ((dpi+dpj)/2.0) - Td) < 0.000001)
                {
                    //Checking interference allowance of most restrictive particle
                    if (dist > (((dpi+dpj)/2.0)-PS) || std::abs(dist-(((dpi+dpj)/2.0)+PS)) < 0.000001)
                    {
                        // Particle is in the safe placement zone (i.e. far enough away from neighbor(i)
                        if (dist < smallestdist)
                        {
                            smallestdist = dist;
                        }
                        neighIndex.push_back(j);
                        neighDist.push_back(dist);
                    }
                    //if else, then the particles are too close, break the loop
                    else
                    {
                        t = 1;                  // Break the for loop, and continue the while loop with t = 1
                        break;
                    }
                }
            }
            
            // Checking to see if the particle is in contact with any other particles
            if (smallestdist > ((dpi+dpj)/2.0))
            {
                // Random value between 1 and 0 inclusive to try the probability of accepting particle placement
                double DispTest = ((double)rand() / (double)RAND_MAX);
                
                // Allowing/Denying particle placement based on probability of agglomerates in dispersion.
                // Probability is based on the "Dispersion Quality" with a perfect (random) dispquality = 1.
                if (DispTest > DispQuality)
                {
                    t = 1; // if the random number is greater than the disp quality, the particle placement is denied, and we break the loop to try again
                    NumDispFail += 1;
                }
            }
            
            // if conditional flag t is still 0, committ the particle coordinates to memory.
            if (t == 0)
            {                
                //Set CBparticle object coordinates
                dat->GetCBParticle(dat->GetParticle(i)->GetTypeID())->SetXYZ(xx, yy, zz, dat->GetXBoundaryLimit(), Td);
                
                //Set maximum prescribable particle diameter for checking edge/face closeness
                int dpm = dat->GetCBMatType(dat->GetParticle(i)->GetType())->GetMaxDiameter();
                
                //Add neighbor indices to particle object
                for (int n = 0; n<neighIndex.size(); n++)
                {
                    dat->GetCBParticle(dat->GetParticle(i)->GetTypeID())->AddNeighbor(neighIndex[n], neighDist[n]);
                    dat->GetCBParticle(dat->GetParticle(neighIndex[n])->GetTypeID())->AddNeighbor(i, neighDist[n]);
                }
                
                //Check closeness to x face boundaries
                if (xx < (LowerX + (dpm + Td)) || std::abs(xx - (LowerX + dpm + Td)) < 0.000001)
                {
                    BoundaryParticles[0].push_back(i);
                }
                else if (xx > (UpperX - (dpm + Td)) || std::abs(xx - (UpperX - (dpm + Td))) < 0.000001)
                {
                    BoundaryParticles[1].push_back(i);
                }
                
                //Check closeness to y face boundaries
                if (yy < (LowerY + (dpm + Td)) || std::abs(yy - (LowerY + (dpm + Td))) < 0.000001)
                {
                    BoundaryParticles[2].push_back(i);
                }
                else if (yy > (UpperY - (dpm+Td)) || std::abs(yy - (UpperY - (dpm + Td))) < 0.000001)
                {
                    BoundaryParticles[3].push_back(i);
                }
                
                //Check closeness to z face boundaries
                if (zz < (LowerZ + (dpm+Td)) || std::abs(zz - (LowerZ + (dpm + Td))) < 0.000001)
                {
                    BoundaryParticles[4].push_back(i);
                }
                else if (zz > (UpperZ - (dpm+Td)) || std::abs(zz - (UpperZ + (dpm + Td))) < 0.000001)
                {
                    BoundaryParticles[5].push_back(i);
                }
            }
            //Clear vectors for next loop
            neighIndex.clear();
            neighDist.clear();
        }
    }
    //Fill corner and edge boundary particle vectors
    CheckCornersAndEdges();
}

void SubRVEdata::CheckCornersAndEdges()
{
    //Filling first edge (x=0, y=0)
    for (int i = 0; i<BoundaryParticles[0].size(); i++)
    {
        if (std::find(BoundaryParticles[2].begin(), BoundaryParticles[2].end(), BoundaryParticles[0][i]) != BoundaryParticles[2].end())
        {
            BoundaryParticles[6].push_back(BoundaryParticles[0][i]);
        }
    }
    //Filling second edge (x=1, y=0)
    for (int i = 0; i<BoundaryParticles[1].size(); i++)
    {
        if (std::find(BoundaryParticles[2].begin(), BoundaryParticles[2].end(), BoundaryParticles[1][i]) != BoundaryParticles[2].end())
        {
            BoundaryParticles[7].push_back(BoundaryParticles[1][i]);
        }
    }
    //Filling third edge (x=1, y=1)
    for (int i = 0; i<BoundaryParticles[1].size(); i++)
    {
        if (std::find(BoundaryParticles[3].begin(), BoundaryParticles[3].end(), BoundaryParticles[1][i]) != BoundaryParticles[3].end())
        {
            BoundaryParticles[8].push_back(BoundaryParticles[1][i]);
        }
    }
    //Filling fourth edge (x=0, y=1)
    for (int i = 0; i<BoundaryParticles[0].size(); i++)
    {
        if (std::find(BoundaryParticles[3].begin(), BoundaryParticles[3].end(), BoundaryParticles[0][i]) != BoundaryParticles[3].end())
        {
            BoundaryParticles[9].push_back(BoundaryParticles[0][i]);
        }
    }
    //Filling the fifth edge (x=0>1, y=0, z=0)
    for (int i = 0; i<BoundaryParticles[2].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[2][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[10].push_back(BoundaryParticles[2][i]);
        }
    }
    //Filling the sixth edge (x=1, y=0>1, z=0)
    for (int i = 0; i<BoundaryParticles[1].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[1][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[11].push_back(BoundaryParticles[1][i]);
        }
    }
    //Filling the seventh edge (x=1>0, y=1, z=0)
    for (int i = 0; i<BoundaryParticles[3].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[3][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[12].push_back(BoundaryParticles[3][i]);
        }
    }
    //Filling the eigth edge (x=0, y=1>0, z=0)
    for (int i = 0; i<BoundaryParticles[0].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[0][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[13].push_back(BoundaryParticles[0][i]);
        }
    }
    //Filling the ninth edge (x=0>1, y=0, z=1)
    for (int i = 0; i<BoundaryParticles[2].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[2][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[14].push_back(BoundaryParticles[2][i]);
        }
    }
    //Filling the tenth edge (x=1, y=0>1, z=1)
    for (int i = 0; i<BoundaryParticles[1].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[1][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[15].push_back(BoundaryParticles[1][i]);
        }
    }
    //Filling the eleventh edge (x=1>0, y=1, z=1)
    for (int i = 0; i<BoundaryParticles[3].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[3][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[16].push_back(BoundaryParticles[3][i]);
        }
    }
    //Filling the twelfth edge (x=0, y=1>0, z=1)
    for (int i = 0; i<BoundaryParticles[0].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[0][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[17].push_back(BoundaryParticles[0][i]);
        }
    }
    //Filling first corner (x=0, y=0, z=0)
    for (int i = 0; i<BoundaryParticles[6].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[6][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[18].push_back(BoundaryParticles[6][i]);
        }
    }
    //Filling second corner (x=1, y=0, z=0)
    for (int i = 0; i<BoundaryParticles[7].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[7][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[19].push_back(BoundaryParticles[7][i]);
        }
    }
    //Filling third corner (x=0, y=1, z=0)
    for (int i = 0; i<BoundaryParticles[9].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[9][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[20].push_back(BoundaryParticles[9][i]);
        }
    }
    //Filling fourth corner (x=1, y=1, z=0)
    for (int i = 0; i<BoundaryParticles[8].size(); i++)
    {
        if (std::find(BoundaryParticles[4].begin(), BoundaryParticles[4].end(), BoundaryParticles[8][i]) != BoundaryParticles[4].end())
        {
            BoundaryParticles[21].push_back(BoundaryParticles[8][i]);
        }
    }
    //Filling fifth corner (x=0, y=0, z=1)
    for (int i = 0; i<BoundaryParticles[6].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[6][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[22].push_back(BoundaryParticles[6][i]);
        }
    }
    //Filling sixth corner (x=1, y=0, z=1)
    for (int i = 0; i<BoundaryParticles[7].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[7][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[23].push_back(BoundaryParticles[7][i]);
        }
    }
    //Filling seventh corner (x=0, y=1, z=1)
    for (int i = 0; i<BoundaryParticles[9].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[9][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[24].push_back(BoundaryParticles[9][i]);
        }
    }
    //Filling eigth corner (x=1, y=1, z=1)
    for (int i = 0; i<BoundaryParticles[8].size(); i++)
    {
        if (std::find(BoundaryParticles[5].begin(), BoundaryParticles[5].end(), BoundaryParticles[8][i]) != BoundaryParticles[5].end())
        {
            BoundaryParticles[25].push_back(BoundaryParticles[8][i]);
        }
    }
}

void SubRVEdata::PlaceCNTParticles(){}

std::vector<std::vector<long int> > &SubRVEdata::GetBoundaryParticles()
{
    return BoundaryParticles;
}

long int SubRVEdata::GetNumDispFail()
{
    return NumDispFail;
}
