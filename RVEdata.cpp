#define _USE_MATH_DEFINES
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <future>
#include <mutex>
#include "Files.h"
#include "PolymerData.h"
#include "CBparticle.h"
#include "CBMaterialData.h"
#include "SubRVEdata.h"
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <random>
#include <Eigen/Core>
#include <Eigen/SparseLU>
#include "RVEdata.h"

RVEdata::RVEdata(){};

RVEdata::~RVEdata(){};

void RVEdata::ReadInputData(std::ifstream &fin, std::ofstream &ferr)
{
    sumInt = 0;
    double NotNeeded;
    std::string lines, ignore;
    
    //Read in basic problem parameters
    fin >> Xbound >> Ybound >> Zbound >> Xsub >> Ysub >> Zsub >> NumCBFillers >> NumCNTFillers >> DispQuality >> d;
    
    //Read in and create CB filler material(s) object(s)
    CB_data.resize(NumCBFillers);
    
    int ID, diamUpper, diamLower, ps;
    float Wcb;
    double rho;
    for (int i = 0; i < NumCBFillers; i++)
    {
        fin >> ID >> Wcb >> rho >> diamLower >> diamUpper >> ps;
        CB_data[i].initialize(ID, rho, diamLower, diamUpper, ps);
        MassFractionOfFillers.push_back(Wcb);
    }
    
    for (int i = 0; i < NumCNTFillers; i++)
    {
        //Read CNT data
    }
    
    long double PolymerBarrierHeight;
    int TunnelingDist;
    float p;
    fin >> rho >> PolymerBarrierHeight >> TunnelingDist >> p;
    p_data.initialize(rho, PolymerBarrierHeight, TunnelingDist, p);
    Td = TunnelingDist;
    
    //Determine quantity of each particle type based on weight percentage
    float Wp = 1;       //polymer at 100% mass fraction
    for (int i = 0; i < MassFractionOfFillers.size(); i++)
    {
        Wp -= MassFractionOfFillers[i]/100.0;      //reduced by mass fraction occupied by filler particles
    }
    long double cells = static_cast<long double>((Xbound - 1))*(Ybound - 1)*(Zbound - 1);
    long double Vrve = cells * (pow(d, 3));     //volume of RVE in nm^3
    
    //Find density of composite using weight fractions https://nptel.ac.in/courses/101106038/mod03lec01.pdf
    double CompRho = Wp/p_data.GetDensity();
    for (int i = 0; i < NumCBFillers; i++)
    {
        CompRho += (MassFractionOfFillers[i]/100.0)/CB_data[i].GetDensity();
    }
    
    //Calculate number of particles of each type and add to CB_data
    for (int i = 0; i<NumCBFillers; i++)
    {
        NumParticles.push_back(round((MassFractionOfFillers[i]/100.0)*Vrve*(1.0/CompRho)/CB_data[i].GetAvgParticleWeight()));
        CB_data[i].SetNumParticles(NumParticles[i]);
    }
    
    for (int i = 0+NumCBFillers; i<NumCBFillers+NumCNTFillers; i++)
    {
        /*
        NumParticles.push_back(round((MassFractionOfFillers[i]/100.0)*Vrve*(1.0/CompRho)/CNT_data[i].GetAvgParticleWeight()));
        CB_data[i].SetNumParticles(NumParticles[i]); */
    }
    
    //Create Sub RVEs
    NumSubRVEs = Xsub*Ysub*Zsub;
    SubRVEs.resize(NumSubRVEs);
    InitializeSubRVEs();
}

void RVEdata::InitializeSubRVEs()
{
    //Defining the volume of the "internal" portion of the sub-RVE's
    double XmajInc = Xbound / Xsub;
    double YmajInc = Ybound / Ysub;
    double ZmajInc = Zbound / Zsub;

    //Check for non integer increments
    if (round(XmajInc) != XmajInc || round(YmajInc) != YmajInc || round(ZmajInc) != ZmajInc)
    {
        std::cout << std::endl << "Non-integer value cannot be assigned to SubRVE Major Increments: Exiting" << std::endl;
        exit(1);
    }
    
    //Determining the Boundary of the internal sub-RVE's and the joining volumes
    int X = 0, Y = 0, Z = 0;
    
    long int SumCBparticles = 0;
    long int SumCNTparticles = 0;
    int CB_count = 0;
    int CNT_count = 0;
    int particle_count = 0;
    std::vector<int> NumSubParticles(NumParticles.size());
    
    Particles NewP;
    CBparticle NewCB;
    CNTparticle NewCNT;
    
    for (int i = 0; i<NumCBFillers; i++)
    {
        NumSubParticles[i] = round(NumParticles[i]/NumSubRVEs);
        SumCBparticles += NumSubParticles[i];
    }
    
    for (int i = NumCBFillers; i < NumCBFillers+NumCNTFillers; i++)
    {
        NumSubParticles[i] = round (NumParticles[i]/NumSubRVEs);
        SumCNTparticles += NumSubParticles[i];
    }
    
    //Create CB particle objects
    for (int i = 0; i<SumCBparticles*NumSubRVEs; i++)
    {
        //Add CB particle with ID count
        CB_particles.push_back(NewCB);
        CB_particles[i].Initialize(CB_count);
        particle_data.push_back(NewP);
        particle_data[particle_count].Initialize(0, CB_particles.size()-1, particle_count);
        CB_count += 1;
        particle_count += 1;
    }
    
    //Create CNT particle objects
    for (int i = 0; i<SumCNTparticles*NumSubRVEs; i++)
    {
        CNT_particles.push_back(NewCNT);
        CNT_particles[i].Initialize(CNT_count);
        particle_data.push_back(NewP);
        particle_data[particle_count].Initialize(1, CNT_particles.size()-1, particle_count);
        CNT_count += 1;
        particle_count += 1;
    }
    
    //Loop to set individual SubRve bounds
    int i = 0;
    int NumTerm = Ysub * Zsub;
    for (int x = 0; x < Xsub; x++)
    {
        for (int z = 0; z < Zsub; z++)
        {
            for (int y = 0; y < Ysub; y++)
            {
                //Initialize SubRVE "i"
                SubRVEs[i].initialize(i, X, X+XmajInc, Y, Y+YmajInc, Z, Z+ZmajInc, NumSubParticles);

                //Increment the subRVE "index" integer and XYZ origin coords.
                i += 1;
                Y += YmajInc;
            }
            Y = 0;
            Z += ZmajInc;
        }
        Z = 0;
        X += XmajInc;
    }
}

int RVEdata::GetNumSubRVEs()
{
    return NumSubRVEs;
}

int RVEdata::GetNumParticles()
{
    return particle_data.size();
}

int RVEdata::GetNumCBFillers()
{
    return NumCBFillers;
}

void RVEdata::FillRVE(int i)
{
    //Place particles in the Sub RVEs
    if (NumCBFillers > 0)
    {
        SubRVEs[i].PlaceCBParticles(i, NumCBFillers, this, Td, d, DispQuality, Xbound);
    }
    if (NumCNTFillers >0)
    {
        SubRVEs[i].PlaceCNTParticles();
    }
}

Particles* RVEdata::GetParticle(int i)
{
    return &particle_data[i];
}

CBparticle* RVEdata::GetCBParticle(int num)
{
    return &CB_particles[num];
}

CNTparticle* RVEdata::GetCNTParticle(int num)
{
    return &CNT_particles[num];
}

CBMaterialData* RVEdata::GetCBMatType(int i)
{
    return &CB_data[i];
}

PolymerData* RVEdata::GetPolymerData()
{
    return &p_data;
}

int RVEdata::GetXBoundaryLimit()
{
    return Xbound;
}

void RVEdata::AssembleRVE()
{   
    //Check interfaces between sub domains for neighboring/interferring particles
    for (int i = 0; i < NumSubRVEs; i++)
    {
        CheckFacesEdgesAndCorners(i);
    }
}

void RVEdata::CheckFacesEdgesAndCorners(int i)
{
    //This function is tricky. Each cubic sub rve has 26 edges, corners, and faces that contain particles in their volue. These particles have the potential to see particles in adjacent sub RVEs. Therefore, we have to "map" each sub RVE by its index, and check the correct Boundary volumes against one another to identify inter-sub-RVE neighbors and interferences.
    //Define vector holding the numbers of the boundary particle vectors the subRve i will need checked against: "Which ones" vector definition: 0 -> subRve i's nth vector, 1 -> adjacent subRve index, 2-> adjacent subRve's jth vector
    std::vector<std::vector<int> > WhichOnes;
    int count = 0;
    
    //The following if statements math the adjacent faces, edges, and corners to be checked against one another. Each set to be checked against is passed into the "WhichOnes" vector (which ones to be checked) and then subsequently particles in those volumes are checked with respect to one another.
    //Check lesser X face if true
    if (i-0.0-(Ysub*Zsub) >= 0)
    {
        WhichOnes.push_back(std::vector<int>(3));
        WhichOnes[count][1] = i-(Ysub*Zsub);
        WhichOnes[count][2] = 1;
        WhichOnes[count][0] = 0;
        count += 1;
        
        //Check lower Z edges if true
        if ((i % (Ysub*Zsub)) >= Ysub)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = i-(Ysub*Zsub)-Ysub;
            WhichOnes[count][2] = 15;
            WhichOnes[count][0] = 13;
            count += 1;
            //Check 1st corner if true
            if (round((i+Ysub+0.0)/Ysub) != (i+Ysub+0.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = i-(Ysub*Zsub)-Ysub-1;
                WhichOnes[count][2] = 25;
                WhichOnes[count][0] = 18;
                count += 1;
            }
            //Check 3rd corner if true
            if (round((i+1.0)/Ysub) != (i+1.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = i-(Ysub*Zsub)-Ysub+1;
                WhichOnes[count][2] = 23;
                WhichOnes[count][0] = 20;
                count += 1;
            }
        }
        
        //Check upper Z edges if true
        if (i % (Ysub*Zsub) < (Ysub*Zsub)-Ysub-0.0)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = (i-(Ysub*Zsub)+Ysub);
            WhichOnes[count][2] = 11;
            WhichOnes[count][0] = 17;
            count += 1;
            //Check 5th Corner if true
            if (round((i+Ysub+0.0)/Ysub) != (i+Ysub+0.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = (i-(Ysub*Zsub)+Ysub-1);
                WhichOnes[count][2] = 21;
                WhichOnes[count][0] = 22;
                count += 1;
            }
            //Check 7th Corner if true
            if (round((i+1.0)/Ysub) != (i+1.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = (i-(Ysub*Zsub)+Ysub+1);
                WhichOnes[count][2] = 19;
                WhichOnes[count][0] = 24;
                count += 1;
            }
        }
    }
    
    //Check greater X face if true
    if (i + 0.0 + (Ysub*Zsub) < NumSubRVEs)
    {
        WhichOnes.push_back(std::vector<int>(3));
        WhichOnes[count][1] = (i+(Ysub*Zsub));
        WhichOnes[count][2] = 0;
        WhichOnes[count][0] = 1;
        count += 1;
        //Check lower edges if true
        if (i % (Ysub*Zsub) >= Ysub)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = (i+(Ysub*Zsub)-Ysub);
            WhichOnes[count][2] = 17;
            WhichOnes[count][0] = 11;
            count += 1;
            //Check 2nd corner if true
            if (round((i+Ysub+ 0.0)/Ysub) != (i+Ysub+ 0.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = (i+(Ysub*Zsub)-Ysub-1);
                WhichOnes[count][2] = 24;
                WhichOnes[count][0] = 19;
                count += 1;
            }
            //Check 4th corner if true
            if (round((i+1.0)/Ysub) != (i+1.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = (i+(Ysub*Zsub)-Ysub+1);
                WhichOnes[count][2] = 22;
                WhichOnes[count][0] = 21;
                count += 1;
            }
        }
        
        //Check upper edges if true
        if (i % (Ysub*Zsub) < (Ysub*Zsub)-Ysub-0.0)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = (i+(Ysub*Zsub)+Ysub);
            WhichOnes[count][2] = 13;
            WhichOnes[count][0] = 15;
            count += 1;
            //Check 6th Corner if true
            if (round((i+Ysub+0.0)/Ysub) != (i+Ysub+0.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = (i+(Ysub*Zsub)+Ysub-1);
                WhichOnes[count][2] = 20;
                WhichOnes[count][0] = 23;
                count += 1;
            }
            //Check 8th Corner if true
            if (round((i+1.0)/Ysub) != (i+1.0)/Ysub)
            {
                WhichOnes.push_back(std::vector<int>(3));
                WhichOnes[count][1] = (i+(Ysub*Zsub)+Ysub+1);
                WhichOnes[count][2] = 18;
                WhichOnes[count][0] = 25;
                count += 1;
            }
        }
    }
    //Check lesser Y face if true
    if (round((i+Ysub+0.0)/Ysub) != (i+Ysub+0.0)/Ysub)
    {
        WhichOnes.push_back(std::vector<int>(3));
        WhichOnes[count][1] = (i-1);
        WhichOnes[count][2] = 3;
        WhichOnes[count][0] = 2;
        count += 1;
        
        //Check lower edges if true
        if (i % (Ysub*Zsub) >= Ysub)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = (i-1-Ysub);
            WhichOnes[count][2] = 16;
            WhichOnes[count][0] = 10;
            count += 1;
        }
        
        //Check upper edges if true
        if (i % (Ysub*Zsub) < ((Ysub*Zsub)-Ysub-0.0))
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = (i-1+Ysub);
            WhichOnes[count][2] = 12;
            WhichOnes[count][0] = 14;
            count += 1;
        }
        
        //Check lesser x-direction edge
        if (i-0.0 - (Ysub*Zsub) >= 0)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = i-(Ysub*Zsub)-1;
            WhichOnes[count][2] = 8;
            WhichOnes[count][0] = 6;
            count += 1;
        }
        
        //Check greater x-direction edge
        if (i + 0.0 + (Ysub*Zsub) < NumSubRVEs)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = i + (Ysub*Zsub) - 1;
            WhichOnes[count][2] = 9;
            WhichOnes[count][0] = 7;
            count += 1;
        }
    }
    
    //Check greater Y face if true
    if (round((i+1.0)/Ysub) != ((i+1.0)/Ysub))
    {
        WhichOnes.push_back(std::vector<int>(3));
        WhichOnes[count][1] = (i+1);
        WhichOnes[count][2] = 2;
        WhichOnes[count][0] = 3;
        count += 1;
        
        //Check lower edges if true
        if (i % (Ysub*Zsub) >= Ysub)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = (i+1-Ysub);
            WhichOnes[count][2] = 14;
            WhichOnes[count][0] = 12;
            count += 1;
        }
        
        //Check upper edges if true
        if (i % (Ysub*Zsub) < (Ysub*Zsub)-Ysub-0.0)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = (i+1+Ysub);
            WhichOnes[count][2] = 10;
            WhichOnes[count][0] = 16;
            count += 1;
        }
        
        //Check lesser x-direction edge
        if ((i-0.0-(Ysub*Zsub)) >= 0)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = i-(Ysub*Zsub)+1;
            WhichOnes[count][2] = 7;
            WhichOnes[count][0] = 9;
            count += 1;
        }
        
        //Check greater x-direction edge
        if (i + 0.0 + (Ysub*Zsub) < NumSubRVEs)
        {
            WhichOnes.push_back(std::vector<int>(3));
            WhichOnes[count][1] = i + (Ysub*Zsub) + 1;
            WhichOnes[count][2] = 6;
            WhichOnes[count][0] = 8;
            count += 1;
        }
    }
    
    //Check lesser Z face if true
    if (i % (Ysub*Zsub) >= Ysub)
    {
        WhichOnes.push_back(std::vector<int>(3));
        WhichOnes[count][1] = i - Ysub;
        WhichOnes[count][2] = 5;
        WhichOnes[count][0] = 4;
        count += 1;
    }
    
    //Check greater Z face if true
    if (i % (Ysub*Zsub) < (Ysub*Zsub)-Ysub-0.0)
    {
        WhichOnes.push_back(std::vector<int>(3));
        WhichOnes[count][1] = i + Ysub;
        WhichOnes[count][2] = 4;
        WhichOnes[count][0] = 5;
        count += 1;
    }

    //Find neighboring particles in the adjacent sub domains
    FindNeighbors(WhichOnes, i);
}

void RVEdata::FindNeighbors(std::vector<std::vector<int> > &WhichOnes, int index)
{
    //This function checks the volumes identified in the check corners and edges function to identify neighbors between the two volumes.
    for (int n = 0; n < WhichOnes.size(); n++)
    {
        //"Which ones" vector definition: 0 -> subRve i's nth vector, 1 -> adjacent subRve index, 2-> adjacent subRve's jth vector
        int P = WhichOnes[n][0];
        int subJ = index;
        int subK = WhichOnes[n][1];
        int Q = WhichOnes[n][2];
        
        for (int j = 0; j < SubRVEs[subJ].GetBoundaryParticles()[P].size(); j++)
        {
            for (int k = 0; k < SubRVEs[subK].GetBoundaryParticles()[Q].size(); k++)
            {
                //boolean that is used to determine whether or not a relationship has already been made between the two particles of interest
                bool AddRelationship;
                double dist;
                //Particle j and k types (CNT, CB, etc)
                int jType, kType;
                //Diameters of particles j and k
                int dpj, dpk;
                //Penetration allowance of j and k
                int jPSoftness, kPSoftness;
                
                //Set values for parameters described above
                jPSoftness = CB_data[CB_particles[particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetTypeID()].GetMatType()].GetPenetrationAllowance();
                kPSoftness = CB_data[CB_particles[particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetTypeID()].GetMatType()].GetPenetrationAllowance();
                dpj = CB_particles[particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetTypeID()].GetDiameter();
                dpk = CB_particles[particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetTypeID()].GetDiameter();
                jType = particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetType();
                kType = particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetType();
                
                //Define the most restrictive penetration allowance
                int PS;
                if (jPSoftness < kPSoftness)
                {
                    PS = jPSoftness;
                }
                else
                {
                    PS = kPSoftness;
                }
                
                //If true, both are CB particles
                if (jType == 0 && kType == 0)
                {
                    dist = CB_particles[particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetTypeID()].FindDistanceToNeighbor(&CB_particles[particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetTypeID()]);
                    //Checking closeness of two particles w.r.t their diameters and tunneling
                    if (dist < ((dpj+dpk)/2.0) + Td || std::abs(dist - ((dpj+dpk)/2.0) - Td) < 0.000001)
                    {
                        if(dist > (((dpj+dpk)/2.0)-PS) || std::abs(dist-(((dpj+dpk)/2.0)-PS)) < 0.000001)
                        {
                            AddRelationship = CheckInVector(CB_particles[particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetTypeID()].GetNeighborsVector(), CB_particles[particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetTypeID()].GetNeighborsVector(), particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetGlobalID(), particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetGlobalID());
                            if (AddRelationship)
                            {
                                CB_particles[particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetTypeID()].AddNeighbor(particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetGlobalID(), dist);
                                CB_particles[particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetTypeID()].AddNeighbor(particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetGlobalID(), dist);
                            }
                        }
                        //if true, particles are interferring at the boundary
                        else if (dist < (((dpj+dpk)/2.0)-PS))
                        {
                            AddRelationship = CheckInVector(CB_particles[particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetTypeID()].GetNeighborsVector(), CB_particles[particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetTypeID()].GetNeighborsVector(), particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetGlobalID(), particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetGlobalID());
                            if (AddRelationship)
                            {
                                CB_particles[particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetTypeID()].AddNeighbor(particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetGlobalID(), dist);
                                CB_particles[particle_data[SubRVEs[subK].GetBoundaryParticles()[Q][k]].GetTypeID()].AddNeighbor(particle_data[SubRVEs[subJ].GetBoundaryParticles()[P][j]].GetGlobalID(), dist);
                                sumInt += 1;
                            }
                        }
                    }
                }
                else if (jType == 1)
                {
                    //Distance between CB and CNT (in progress)
                }
                else if (kType == 1)
                {
                    //Distance between CB and CNT (in progress)
                }
            }
        }
    }
}

void RVEdata::SolveForConductivity()
{
    //Vectors passed to the find network algorithm
    long int N = particle_data.size();
    
    std::mutex *PosLocks;
    std::mutex *GndLocks;
    
    PosLocks = new std::mutex[N];
    GndLocks = new std::mutex[N];
    
    posbank.resize(N, std::vector<bool> (2));
    gndbank.resize(N, std::vector<bool> (2));
    network.resize(N);
    
    //Initialize vector lock objects
    for (long int i = 0; i<N; i++)
    {
        posbank[i][0] = false;
        gndbank[i][0] = false;
        posbank[i][1] = false;
        gndbank[i][1] = false;
    }
    
    //Find the networked particles
    NetExists = FindNetwork(N, &PosLocks, &GndLocks);
    
    if (NetExists == false)
    {
        return;
    }
    
    //If arrived at this point, solve for conducitiviy of network particles
    SolveForNodalVoltages();

    RVEresistance = VIN/sumCurrent;
    RVEconductivity = (Xbound*d)/(Ybound*Zbound*d*d*RVEresistance);
}

bool RVEdata::FindNetwork(long int N, std::mutex *PosLocks[], std::mutex *GndLocks[])
{
    //Fill initial posbank and gndbank first rows
    for (long int i = 0; i<particle_data.size(); i++)
    {
        if (particle_data[i].GetType() < NumCBFillers)
        {
            //It's a cb particle
            posbank[i][0] = CB_particles[particle_data[i].GetTypeID()].IsPositive();
            gndbank[i][0] = CB_particles[particle_data[i].GetTypeID()].IsNegative();
        }
        else
        {
            //It's a cnt particle
        }
    }
    
    //Count num of processesors available
    int NumProcessors = std::thread::hardware_concurrency();

    //Half of available processors do positive search, half do negative search
    std::vector<std::future<void> > BranchSearchSplit;
    
    BranchSearchSplit.push_back(std::async(std::launch::async, &RVEdata::PositiveBranch, this, NumProcessors/2, PosLocks));
    BranchSearchSplit.push_back(std::async(std::launch::async, &RVEdata::NegativeBranch, this, NumProcessors/2, GndLocks));
    BranchSearchSplit[0].wait();
    BranchSearchSplit[1].wait();
    
    //Check for Network
    int NumNetwork = 0;
    
    for (unsigned long int i = 0; i<N; i++)
    {
        if (posbank[i][1] == true && gndbank[i][1] == true)
        {
            network[i] = 1;
            NumNetwork += 1;
        }
        else
        {
            network[i] = 0;
        }
    }
    
    if (NumNetwork > 0)
    {
        return(true);
    }
    else
    {
        return(false);
    }
}

void RVEdata::PositiveBranch(int NumProcessors, std::mutex *PosLocks[])
{
    //initializing variables that help identify when all pos network has been found
    int NumPos2 = 0;
    int NumPos = 1;
    
    // This while loop is dependent on the two positive network identifiers. In each iteration
    // the code will find positive network particles identified by a 1 in the first row of the
    // posbank matrix. If the particle has a "0" in the second row and a 1 in the first row,
    // this signals that the particle has not yet been processed. The processing involves searching
    // the particle for its neighbors, and changing those neighbor indexes in the 1st row of the
    // posbank to 1's. This process is repeated until there are no more changes in the number of
    // 1's in the posbank from the beginning of the loop to the end of the loop.
    
    while (NumPos != NumPos2)
    {
        // Finding the total number of particles in the positive network at the loop start
        NumPos = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (posbank[w][1] == true)
            {
                NumPos += 1;
            }
        }
        
        //Split the search party into futures
        std::vector<std::future<void> > SearchersPos;
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersPos.push_back(std::async(std::launch::async, &RVEdata::PositiveSearch, this, PosLocks));
        }
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersPos[i].wait();
        }
        
        // Finding the total number of particles in the pos network at the end of the loop
        NumPos2 = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (posbank[w][1] == true)
            {
                NumPos2 += 1;
            }
        }
    }
    std::cout << "NumPos = " << NumPos2 << std::endl;
}

//Algorithm that searches the posbank for particles that need to be searched and have their neighbors added to posbank
void RVEdata::PositiveSearch(std::mutex *PosLocks[])
{
    //Neighbors vector placeholder
    std::vector<int> neighbors;
    unsigned long int NumNeigh = 0;
    
    for (unsigned long int i = 0; i<particle_data.size(); i++)
    {
        if (PosLocks[i]->try_lock() == true)
        {
            if (posbank[i][0] == true && posbank[i][1] == false)
            {
                neighbors = GetNeighbors(i);
                NumNeigh = neighbors.size();
                if (NumNeigh != 0)
                {
                    for (unsigned int b = 0; b<NumNeigh; b++)
                    {
                        bool Success = false;
                        int FailCount = 0;
                        while (Success == false)
                        {
                            if (PosLocks[neighbors[b]]->try_lock() == true)
                            {
                                posbank[neighbors[b]][0] = true;
                                PosLocks[neighbors[b]]->unlock();
                                Success = true;
                            }
                            else
                            {
                                FailCount += 1;
                            }
                            
                            if (FailCount > 5)
                            {
                                PosLocks[i]->unlock();
                                goto TryLaterNeg;
                            }
                        }
                    }
                }
                posbank[i][1] = 1;
            }
            PosLocks[i]->unlock();
        TryLaterNeg:
            continue;
        }
        else
        {
            continue;
        }
    }
}

void RVEdata::NegativeBranch(int NumProcessors, std::mutex *GndLocks[])
{
    // Finding the total number of particles direclty connected to positive electrode
    int NumGnd = 1;
    int NumGnd2 = 0;
    // This while loop is dependent on the two negative network identifiers. In each iteration
    // the code will find negative network particles identified by a 1 in the first row of the
    // gndbank matrix. If the particle has a "0" in the second row and a 1 in the first row,
    // this signals that the particle has not yet been processed. The processing involves searching
    // the particle for its neighbors, and changing those neighbor indexes in the 1st row of the
    // gndbank to 1's. This process is repeated until there are no more changes in the number of
    // 1's in the gndbank from the beginning of the loop to the end of the loop.

    while (NumGnd != NumGnd2)
    {
        // Finding the total number of particles in the ground network at the loop start
        NumGnd = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (gndbank[w][1] == true)
            {
                NumGnd += 1;
            }
        }
        
        //Split the search party into futures
        std::vector<std::future<void> > SearchersNeg;
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersNeg.push_back(std::async(std::launch::async, &RVEdata::NegativeSearch, this, GndLocks));
        }
        
        for (int i = 0; i<NumProcessors; i++)
        {
            SearchersNeg[i].wait();
        }
        
        // Finding the total number of particles in the ground network at the loop end
        NumGnd2 = 0;
        for (unsigned long int w = 0; w<particle_data.size(); w++)
        {
            if (gndbank[w][1] == true)
            {
                NumGnd2 += 1;
            }
        }
    }
    std::cout << "NumGnd = " << NumGnd2 << std::endl;
}

//Algorithm that searches the gndbank for particles that need to be searched and have their neighbors added to gndbank
void RVEdata::NegativeSearch(std::mutex *GndLocks[])
{
    //Neighbors vector placeholder
    std::vector<int> neighbors;
    unsigned long int NumNeigh = 0;

    for (unsigned long int i = 0; i<particle_data.size(); i++)
    {
        if (GndLocks[i]->try_lock() == true)
        {
            if (gndbank[i][0] == true && gndbank[i][1] == false)
            {
                neighbors = GetNeighbors(i);
                NumNeigh = neighbors.size();
                if (NumNeigh != 0)
                {
                    for (unsigned int b = 0; b<NumNeigh; b++)
                    {
                        bool Success = false;
                        int FailCount = 0;
                        while (Success == false)
                        {
                            if (GndLocks[neighbors[b]]->try_lock() == true)
                            {
                                gndbank[neighbors[b]][0] = true;
                                GndLocks[neighbors[b]]->unlock();
                                Success = true;
                            }
                            else
                            {
                                FailCount += 1;
                            }
                            
                            if (FailCount > 5)
                            {
                                GndLocks[i]->unlock();
                                goto TryLaterNeg;
                            }
                        }
                    }
                }
                gndbank[i][1] = 1;
            }
            GndLocks[i]->unlock();
        TryLaterNeg:
            continue;
        }
        else
        {
            continue;
        }
    }
}

void RVEdata::OutputResults(std::ofstream &fout)
{
    long int TotalDispFails = 0;
    for (int i = 0; i < NumSubRVEs; i++)
    {
        TotalDispFails += SubRVEs[i].GetNumDispFail();
    }
    
    fout << RVEresistance << "," << RVEconductivity << "," << (static_cast<double>(sumInt)*2.0/static_cast<double>(particle_data.size())) << "," << TotalDispFails << ",";
}

bool RVEdata::CheckInVector(std::vector<int> &ParticleN, std::vector<int> &ParticleO, int ParticleIndexN, int ParticleIndexO)
{
    if (std::find(ParticleN.begin(), ParticleN.end(), ParticleIndexO) != ParticleN.end())
    {
        return (false);
    }
    else if (std::find(ParticleO.begin(), ParticleO.end(), ParticleIndexN) != ParticleO.end())
    {
        return (false);
    }
    else
    {
        return (true);
    }
}

std::vector<int> &RVEdata::GetNeighbors(int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        return CB_particles[particle_data[i].GetTypeID()].GetNeighborsVector();
    }
    else
    {
        //Do nothing
        std::cout << "ERROR: CNT particles don't exist yet" << std::endl;
        exit(1);
    }
}

std::vector<double> &RVEdata::GetNeighborsDist(int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        return CB_particles[particle_data[i].GetTypeID()].GetNeighborsDistVector();
    }
    else
    {
        //Do nothing
        std::cout << "ERROR: CNT particles don't exist yet" << std::endl;
        exit(1);
    }
}

void RVEdata::SolveForNodalVoltages()
{
    int VIN = 5;
    //Count num particles in network and assign row indices to network particles
    std::vector<int> Indices;
    std::vector<int> NeighborIndices;
    int SumNet = 0;
    for (int i = 0; i < network.size(); i++)
    {
        if (network[i] == true)
        {
            Indices.push_back(i);
            NeighborIndices.push_back(SumNet);
            SumNet += 1;
        }
        else
        {
            NeighborIndices.push_back(-1);
        }
    }
    
    double Rc = 1000;          //Contact resistance between carbon allotropes in ohms
    
    Eigen::SparseMatrix<double, Eigen::ColMajor> g (SumNet, SumNet);
    g.reserve(10);
    Eigen::VectorXd curr = Eigen::VectorXd::Zero(SumNet);
    
    //Find the particle distribution between cores
    int NumP = SumNet / NumSubRVEs;
    int Left = SumNet % NumSubRVEs;
    
    //Place particles in the SubRVEs using threads
    std::vector<std::future<std::vector<Eigen::Triplet<double> > > > Waiters;
    typedef Eigen::Triplet<double> T;
    std::vector<T> GtripletList;
    
    std::vector<std::vector<int> > IndicesG;
    std::vector<std::vector<double> > ValuesG;
    std::vector<double> ValuesCurr;
    IndicesG.resize(SumNet);
    ValuesG.resize(SumNet);
    ValuesCurr.resize(SumNet);
    
    int start = 0;
    int finish;
    
    auto FillStart = std::chrono::high_resolution_clock::now();
    std::cout << "Filling Cond " << std::flush;
    
    //Loop to fill G and Curr using multiple threads
    for (int i = 0; i<NumSubRVEs; i++)
    {
        if (i == NumSubRVEs - 1)
        {
            finish = SumNet;
        }
        else
        {
            finish = start + NumP;
        }
        Waiters.push_back(std::async(std::launch::async, &RVEdata::FillCurrAndG, this, std::ref(ValuesCurr), std::ref(Indices), std::ref(NeighborIndices), Rc, VIN, start, finish));
        start += NumP;
    }
    
    std::vector<Eigen::Triplet<double> > Gtemp;
    
   for (int i = 0; i<Waiters.size(); i++)
    {
        Waiters[i].wait();
        Gtemp = Waiters[i].get();
        GtripletList.insert(GtripletList.end(), Gtemp.begin(), Gtemp.end());
        Gtemp.clear();
    }
    
    //Fill G and Curr Eigen objects using triplet list
    for (int i = 0; i<SumNet; i++)
    {
        //Current vector portion
        if (ValuesCurr[i] > 0)
        {
            curr(i) = ValuesCurr[i];
        }
    }
    
    g.setFromTriplets(GtripletList.begin(), GtripletList.end());
    
    auto FillFinish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> FilledElapsed = FillFinish - FillStart;
    std::cout << FilledElapsed.count() << "..Solving Linear Alg " << std::flush;

    //Perform linear algebra to obtain nodal voltages
    g.finalize();
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > solver;
    g.makeCompressed();
    solver.analyzePattern(g);
    solver.compute(g);
    Eigen::VectorXd V = Eigen::VectorXd::Zero(SumNet,1);
    solver.factorize(g);
	
    if (solver.info() != 0)
    {
        std::cout << std::endl << "Matrix Factorization Failed!" << std::endl<< solver.lastErrorMessage() << std::endl;
        exit (1);
    }
    V = solver.solve(curr);
    
    CheckNeighbors();
    auto Solved = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> SolvedFinal = Solved - FillFinish;
    std::cout << SolvedFinal.count() << ".." << std::flush;
    
    //Fill Voltage vector in RVEdata class with nodal voltages found
    Nodal_Voltages.resize(particle_data.size());
    int ph = 0;
    for (int i = 0; i < particle_data.size(); i++)
    {
        if (network[i] == false)
        {
            Nodal_Voltages[i] = 0;
            continue;
        }
        Nodal_Voltages[i] = V.coeffRef(ph);
        ph += 1;
        
        if (Nodal_Voltages[i]-VIN > 0.001)
        {
            std::cout << "i: " << i << " V: " << Nodal_Voltages[i] << std::endl;
        }
    }
    
    //Sum the current flowing into the positively connected particles and calculate resistance/conductivity
    for (int i = 0; i < SumNet; i++)
    {
        if (particle_data[Indices[i]].GetPosTerminalResistance() != 0)
        {
            sumCurrent += (VIN - Nodal_Voltages[Indices[i]])/particle_data[Indices[i]].GetPosTerminalResistance();
        }
    }
}

std::vector<Eigen::Triplet<double> > RVEdata::FillCurrAndG(std::vector<double> &ValuesCurr, std::vector<int> &Indices, std::vector<int> &NeighborIndices, double Rc, int VIN, int start, int finish)
{
    //Neighbor vector placeholder and distance placeholder and tunneling threshold
    std::vector<int> neighbors;
    std::vector<double> neighborsDist;
    std::vector<Eigen::Triplet<double> > G;
    double TunnelingThresh;
    double TunnelingResistance;
    double ValuesG;
    
    for (int i = start; i<finish; i++)
    {        
        //Fill Conductivity matrix using data from neighbor vector
        neighbors = GetNeighbors(Indices[i]);
        neighborsDist = GetNeighborsDist(Indices[i]);
        for (int j = 0; j < neighbors.size(); j++)
        {
            TunnelingThresh = CalcTunnelingThreshold(Indices[i],neighbors[j]);
            if (neighborsDist[j] < TunnelingThresh || std::abs(neighborsDist[j] - TunnelingThresh) < 0.0001)
            {
                //Particles are in CONTACT
                ValuesG = -1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,NeighborIndices[neighbors[j]],ValuesG));
                ValuesG = 1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
            else
            {
                long double At = FindTunnelingArea(Indices[i],neighbors[j]);
                TunnelingResistance = CalcTunnelingResistance(neighborsDist[j], TunnelingThresh, d, At);
                if (TunnelingResistance < Rc)
                {
                    TunnelingResistance = Rc;
                }
                ValuesG = -1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,NeighborIndices[neighbors[j]],ValuesG));
                ValuesG = 1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
        }
        //Check particle i's relation to positive and negative terminals
        if (CheckPosNeg(Indices[i],false) == true)
        {
            //Find distance to positive terminal and place conductivity in g matrix and fill current vector with corresponding value
            double dist = FindDistanceToTerminal(false,Indices[i]);
            TunnelingThresh = CalcTunnelingThresholdToTerminal(Indices[i]);
            if (dist < TunnelingThresh || std::abs(dist - TunnelingThresh) < 0.0001)
            {
                //Particles are in CONTACT
                particle_data[Indices[i]].SetPosTerminalResistance(Rc);
                ValuesCurr[i] = (VIN/Rc);
                ValuesG = 1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
            else
            {
                //Particles are tunneling
                long double At = FindTunnelingArea(Indices[i],-1);
                TunnelingResistance = CalcTunnelingResistance(dist, TunnelingThresh, d, At);
                if (TunnelingResistance < Rc)
                {
                    TunnelingResistance = Rc;
                }
                
                particle_data[Indices[i]].SetPosTerminalResistance(TunnelingResistance);
                ValuesCurr[i] = (VIN/TunnelingResistance);
                ValuesG = 1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
        }
        if (CheckPosNeg(Indices[i],true) == true)
        {
            //Find distance to negative terminal and place conductivity in g matrix
            double dist = FindDistanceToTerminal(true,Indices[i]);
            TunnelingThresh = CalcTunnelingThresholdToTerminal(Indices[i]);
            if (dist < TunnelingThresh || std::abs(dist - TunnelingThresh) < 0.0001)
            {
                //Particles are in CONTACT
                ValuesG = 1.0/Rc;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
            else
            {
                //particles are tunneling
                long double At = FindTunnelingArea(Indices[i],-1);
                TunnelingResistance = CalcTunnelingResistance(dist, TunnelingThresh, d, At);
                if (TunnelingResistance < Rc)
                {
                    TunnelingResistance = Rc;
                }
                ValuesG = 1.0/TunnelingResistance;
                G.push_back(Eigen::Triplet<double>(i,i,ValuesG));
            }
        }
    }
    return G;
}

double RVEdata::CalcTunnelingThreshold(int i, int j)
{
    int iType, jType;
    iType = particle_data[i].GetType();
    jType = particle_data[j].GetType();

    if (iType < NumCBFillers && jType < NumCBFillers)
    {
        //Both are CB particles
        int dpi, dpj;
        dpi = CB_particles[particle_data[i].GetTypeID()].GetDiameter();
        dpj = CB_particles[particle_data[j].GetTypeID()].GetDiameter();
        return ((dpi+dpj)/2.0);
    }
    else if (iType > NumCBFillers && jType > NumCBFillers)
    {
        //Both i and j are CNT particles
        return 0;
    }
    else if (iType > NumCBFillers && jType < NumCBFillers)
    {
        //i is a cnt and j is a cb
        return 0;
    }
    else if (iType < NumCBFillers && jType > NumCBFillers)
    {
        //j is a cnt and i is a cb
        return 0;
    }
    else
    {
        return 0;
    }
}

double RVEdata::CalcTunnelingThresholdToTerminal(int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        //i is a cb particle
        return (CB_particles[particle_data[i].GetTypeID()].GetDiameter()/2.0);
    }
    else
    {
        return 0;
    }
}

double RVEdata::CalcTunnelingResistance(double dist, double TunnelingThresh, double u, long double At)
{
    //Tunneling Conductivity Parameter definition
    long double h = 6.626*pow(10,-34);           //Planck's constant (m^2kg/s)
    long double e = 1.602176634*pow(10,-19);     //Elementary charge in units of Coulombs
    long double tau = 0.5*1.60218*pow(10,-19);   //Average barrier height of Epoxy
    long double Em = 9.10938356*pow(10,-31);     //electron mass
    double d;                                   //tunneling distance
    
    d = (dist - TunnelingThresh)*u;
    
    if (d < 0)
    {
        std::cout << "Can't have negative tunneling distance.." << std::endl;
        exit(1);
    }
    
    long double RtunnNumerator = exp(4.0*M_PI*d*sqrt(2.0 * Em*tau) / h)*(h*h*d);
    long double RtunnDenom = (At*e*e*sqrt(2.0 * Em*tau));
    
    double TR = RtunnNumerator/RtunnDenom;
    
    if (dist - TunnelingThresh > 2)
    {
        std::cout << "Why? How? Dist = " << dist << " TT = " << TunnelingThresh << std::endl;
    }
    
    if (TR == 0)
    {
        std::cout << "Tunn R can't be 0: Rnumerator = " << RtunnNumerator << " dist= " << dist << " thresh= " << TunnelingThresh << " d= " << d << " u = " << u << " 2.0*EM*tau= " << 2.0*Em*tau << std::endl;
    }
    else if (TR < 0)
    {
        std::cout << "Tunn R can't be less than 0: Rnumer = " << RtunnNumerator << " Rdenom= " << RtunnDenom << std::endl;
    }
    if (At == 0)
    {
        std::cout << " At: " << At;
    }
    return TR;
}

double  RVEdata::FindDistanceToTerminal(bool terminal, int i)
{
    if (particle_data[i].GetType() < NumCBFillers)
    {
        //i is a cb particle
        if (terminal == false)
        {
            //positive terminal
            return CB_particles[particle_data[i].GetTypeID()].GetXCoordinate();
        }
        else
        {
            //ground terminal
            return Xbound - CB_particles[particle_data[i].GetTypeID()].GetXCoordinate();
        }
    }
    else
    {
        return 0;
    }
}

long double RVEdata::FindTunnelingArea(int i, int j)
{
    if (j < 0)
    {
        //Particle to terminal tunneling
        if (particle_data[i].GetType() <= NumCBFillers)
        {
            //CB particle to terminal
            double dp = CB_particles[particle_data[i].GetTypeID()].GetDiameter()*pow(10,-9);
            return M_PI*pow(dp/2.0,2);           //cross sectional area of tunneling path (m^2)
        }
        else
        {
            //CNT particle to terminal
            return 0;
        }
    }
    if (j >= 0)
    {
        //Particle to particle tunneling
        if (particle_data[i].GetType() <= NumCBFillers && particle_data[j].GetType() <= NumCBFillers)
        {
            //CB to CB tunneling
            double dpi = CB_particles[particle_data[i].GetTypeID()].GetDiameter()*pow(10,-9);
            double dpj = CB_particles[particle_data[j].GetTypeID()].GetDiameter()*pow(10,-9);
            return M_PI*pow((dpi+dpj)/4.0,2);           //cross sectional area of tunneling path (m^2)
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

bool RVEdata::CheckPosNeg(int i, bool j)
{
    int Type = particle_data[i].GetType();
    
    if (j == true)
    {
        //Check to Negative terminal
        if (Type < NumCBFillers)
        {
            //CB particle
            return CB_particles[particle_data[i].GetTypeID()].IsNegative();
        }
        else
        {
            //Not CB particle
            return 0;
        }
    }
    if (j == false)
    {
        //Check to positive terminal
        if (Type < NumCBFillers)
        {
            //CB particle
            return CB_particles[particle_data[i].GetTypeID()].IsPositive();
        }
        else
        {
            //Not CB particle
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

void RVEdata::CheckNeighbors()
{
    std::string dir = getcwd(NULL, 0);
    std::ofstream MFErr;
    MFErr.open(dir + "/CheckNeighbors.csv");
    for (int i = 0; i < particle_data.size(); i++)
    {
        std::vector<int> neigh = CB_particles[i].GetNeighborsVector();
        MFErr << CB_particles[i].GetXCoordinate() << "," << CB_particles[i].GetYCoordinate() << "," << CB_particles[i].GetZCoordinate() << ",";
        if (network[i] == true)
        {
            MFErr << "1,Neighbors,";
        }
        else
        {
            MFErr << "0,Neighbors,";
        }
        for (int j = 0; j < neigh.size(); j++)
        {
            MFErr << neigh[j] << ",";
        }
        MFErr << std::endl;
    }
    MFErr.close();
}
