#include "CNTparticle.h"

CNTparticle::CNTparticle(){}

CNTparticle::~CNTparticle(){}

void CNTparticle::Initialize(int i)
{
    ID = i;
}

int CNTparticle::GetID()
{
    return ID;
}
