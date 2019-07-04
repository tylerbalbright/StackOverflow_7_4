//This class holds the particle information such as particle type, and subsequent particle ID based on that type specifier.

#pragma once

class Particles
{
private:
    int Type;
    int TypeID;
    int GlobalID;
    double PosTerminalResistance;
    
public:
    //Constructor and Destructor
    Particles();
    ~Particles();
    
    void Initialize(int t, int tID, int gID);
    int GetType();
    int GetTypeID();
    int GetGlobalID();
    void SetPosTerminalResistance(double r);
    double GetPosTerminalResistance();
};
