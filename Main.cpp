#include <iostream>
#include <fstream>
#include <vector>
#include <future>
#include <chrono>
#include "CBMaterialData.h"
#include "PolymerData.h"
#include "Files.h"
#include "CBparticle.h"
#include "RVEdata.h"
#include "SubRVEdata.h"

int main()
{
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    
    std::ifstream fin;      //stream for problem input
    std::ofstream fout;      //stream for problem output
    std::ofstream ferr;      //stream for error output
    std::ofstream fplot;     //stream for plotting RVE
    
    RVEdata dat;        //object to store problem data
    Files o_files;      //object to manipulate files
    
    //Open all file streams
    o_files.OpenFiles(fin, fout, fplot, ferr);
    
    //Read in data for problem
    dat.ReadInputData(fin, ferr);
    fin.close();
    std::cout << "Placing " << std::flush;
    //Place particles in the SubRVEs using threads
    std::vector<std::future<void> > Waiters;
    for (int i = 0; i < dat.GetNumSubRVEs(); i++)
    {
        Waiters.push_back(std::async(std::launch::async,&RVEdata::FillRVE, &dat, i));
    }
    for (int i = 0; i < Waiters.size(); i++)
    {
        Waiters[i].wait();
    }
    auto placed = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> placedElapsed = placed - start;
    std::cout << placedElapsed.count() << ".." << std::flush;

    std::ofstream XYZ;
    XYZ.open("XYZ.csv");
    XYZ << dat.GetNumParticles() << std::endl << std::endl;
    for (int k = 0; k<dat.GetNumParticles(); k++)
    {
        XYZ << dat.GetCBParticle(k)->GetXCoordinate() << "," << dat.GetCBParticle(k)->GetYCoordinate() << "," << dat.GetCBParticle(k)->GetZCoordinate() << std::endl;
    }
    XYZ.close(); 
    
    std::cout << "Assembling " << std::flush;
    //Assemble full RVE
    dat.AssembleRVE();
    
    auto Assembled = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> AssembledElapsed = Assembled - placed;
    std::cout << AssembledElapsed.count() << ".." << std::flush;
    
    //Solve for RVE conductivity
    dat.SolveForConductivity();
    
    //Output result files
    dat.OutputResults(fout);
    
    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    
    fout << elapsed.count() << std::endl;
    
    //Clos all files
    o_files.CloseFiles(fout, fplot, ferr);
    std::cout << "Done.." << std::endl;
}
