/* The files class handles the file storage functions of input and output data. */
#include <string>

#pragma once

class Files
{
private:
    std::string input;
    std::string error;
    std::string output;
    std::string OvitoPlot;
    
public:
    //Constructor and Destructor
    Files();
    ~Files();
    
    //Functions to open and close files
    void OpenFiles(std::ifstream &fin, std::ofstream &fout, std::ofstream &fplot, std::ofstream &ferr);
    void CloseFiles(std::ofstream &fout, std::ofstream &fplot, std::ofstream &ferr);
};
