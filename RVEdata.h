/* The RVE data class holds all of the data and major functions required to store the runtime data. */
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <mutex>
#include "Particles.h"
#include "CBparticle.h"
#include "CNTparticle.h"
#include "CBMaterialData.h"
#include "PolymerData.h"
#pragma once

class SubRVEdata;
class RVEdata
{
private:
    int Xbound, Ybound, Zbound;   //RVE size definition
    int NumCBFillers;             //Number of CB filler types in RVE
    int NumCNTFillers;           //Number of CNT filler types in RVE
    int Xsub, Ysub, Zsub;       //Defining parallelization volume counts in the x, y, and z axes.
    int NumSubRVEs;             //Defining total numbe of sub RVEs in the parallelized RVE
    int Td;                     //Tunneling distance in nanometers
    double d;                  //unit size in nanometers
    float DispQuality;    //Float value between zero and one
    int sumInt = 0;     //Total interferring particles in dataset
    double sumCurrent = 0;  //Summing the current flowing into the positively connected particles
    std::vector<bool> network;
    std::vector<std::vector<bool> > posbank;
    std::vector<std::vector<bool> > gndbank;
    double RVEresistance;   //Resistance across the RVE network
    double RVEconductivity; //Conductivity of the RVE network
    double VIN = 5;     //Voltage differential applied across the RVE
    bool NetExists;     //Boolean to flag if a network exists in the RVE
    std::vector<float> MassFractionOfFillers;     //Define volume fractions of filler(s)
    std::vector<SubRVEdata> SubRVEs;        //Array of SubRVE data storage objects
    std::vector<CBMaterialData> CB_data;     //Array of CBMaterialData storage objects
    std::vector<CBparticle> CB_particles;  //2-D Array of CBparticle types storage objects
    std::vector<CNTparticle> CNT_particles;  //2-D Array of CNT particle types storage objects
    PolymerData p_data;         ///Object storing info about polymer
    std::vector<long int> NumParticles;  //vector storing number of each type of particle in RVE
    std::vector<std::vector<int> > BoundaryParticles;  //Vector storing objects that have interfernce potential
    std::vector<Particles> particle_data;       //Vector of particle objects holding details about the particle types and their local/global IDs
    std::vector<double> Nodal_Voltages;
    std::vector<std::vector<double> > ConductivityMatrix;

public:
    //Constructor and Destructor
    RVEdata();
    ~RVEdata();
    
    //Function to read input file
    void ReadInputData(std::ifstream &fin, std::ofstream &ferr);
    
    //Initialize the sub domains to distribute the computing load
    void InitializeSubRVEs();
    
    //Return number of sub domains in problem
    int GetNumSubRVEs();
    
    //Return positive-ground distance (limit)
    int GetXBoundaryLimit();
    
    //Populate the RVE with the prescribed particles defined in the input file
    void FillRVE(int i);
    
    //Return number of particles
    int GetNumParticles();
    
    //Return number of CB fillers
    int GetNumCBFillers();
    
    //Return point to particle i
    Particles* GetParticle(int i);

    //Return pointer to CB particle #i
    CBparticle* GetCBParticle(int num);
    
    //Return pointer to CNT particle #i
    CNTparticle* GetCNTParticle(int num);
    
    //Return pointer to material type object
    CBMaterialData* GetCBMatType(int i);
    
    //Return pointer to polymer data object
    PolymerData* GetPolymerData();
    
    //Function to assemble sub RVEs into main RVE
    void AssembleRVE();
    
    //Find particle neighbors in coincident RVEs
    void FindNeighbors(std::vector<std::vector<int> > &WhichOnes, int i);
    
    //Functions to check coincident volumes of sub domains for neighboring particles
    void CheckFacesEdgesAndCorners(int i);
    
    //Check vector for value N in vector Vo and value O in vector Vn
    bool CheckInVector(std::vector<int> &ParticleN, std::vector<int> &ParticleO, int ParticleIndexN, int ParticleIndexO);
    
    //Function to calculate the conductivity of the RVE
    void SolveForConductivity();
    
    //Function to solve for the nodal voltages of the networked RVE particles
    void SolveForNodalVoltages();
    
    //Function to thread the filling of the conductivity and current matrices
    std::vector<Eigen::Triplet<double> > FillCurrAndG(std::vector<double> &ValuesCurr, std::vector<int> &Indices, std::vector<int> &NeighborIndices, double Rc, int VIN, int start, int finish);

    //Find the network via positive -> negative algorithm
    bool FindNetwork(long int N, std::mutex *PosLocks[], std::mutex *GndLocks[]);
    
    //Branch searching functions for parrallel find net algorithm
    void NegativeBranch(int NumProcessors, std::mutex *GndLocks[]);
    void PositiveBranch(int NumProcessors, std::mutex *PosLocks[]);
    void PositiveSearch(std::mutex *PosLocks[]);
    void NegativeSearch(std::mutex *GndLocks[]);

    //Return instance of neighbors vector for particle i
    std::vector<int> &GetNeighbors(int i);
    
    //Return instance of neighbors distance vector for particle i
    std::vector<double> &GetNeighborsDist(int i);

    //Return tunneling threshold between two particles
    double CalcTunnelingThreshold(int i, int j);
    
    //Return tunneling threshold between particle and terminal
    double CalcTunnelingThresholdToTerminal(int i);
    
    //Return cross-sectional area of tunneling path
    long double FindTunnelingArea(int i, int j);

    //Return tunneling resistance between two particles
    double CalcTunnelingResistance(double dist, double TunnelingThresh, double u, long double At);
    
    //Return distance from particle to terminal (pos or gnd defined by boolean)
    double  FindDistanceToTerminal(bool terminal, int i);
    
    //Check if particle is close to pos or neg terminal
    bool CheckPosNeg(int i, bool j);
    
    //Error reporting
    void CheckNeighbors();
    
    void OutputResults(std::ofstream &fout);
};
