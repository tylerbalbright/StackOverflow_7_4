#include <fstream>
#include <string>
#include <iostream>
#include <unistd.h>
#include "Files.h"

Files::Files()
{
    std::string cwd;
    cwd = getcwd(NULL,0);
    
    input = (cwd + "/Input_Parameters.txt");
    output = (cwd + "/Output.csv");
    error = (cwd + "/error.csv");
    OvitoPlot = (cwd + "/OvitoPlot.csv");
}

Files::~Files(){}

void Files::OpenFiles(std::ifstream &fin, std::ofstream &fout, std::ofstream &fplot, std::ofstream &ferr)
{
    fin.open(input);
    fout.open(output, std::ios::app);
    //fplot.open(OvitoPlot);
    ferr.open(error);
}

void Files::CloseFiles(std::ofstream &fout, std::ofstream &fplot, std::ofstream &ferr)
{
    fout.close();
    //fplot.close();
    ferr.close();
}

