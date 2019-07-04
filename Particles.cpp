#include "Particles.h"

Particles::Particles()
{
    Type = -1;
    TypeID = -1;
    GlobalID = -1;
    PosTerminalResistance = 0;
}
Particles::~Particles(){}

void Particles::Initialize(int t, int tID, int gID)
{
    Type = t;
    TypeID = tID;
    GlobalID = gID;
    PosTerminalResistance = 0;
}

int Particles::GetType()
{
    return Type;
}

int Particles::GetTypeID()
{
    return TypeID;
}

int Particles::GetGlobalID()
{
    return GlobalID;
}

double Particles::GetPosTerminalResistance()
{
    return PosTerminalResistance;
}

void Particles::SetPosTerminalResistance(double r)
{
    PosTerminalResistance = r;
}
