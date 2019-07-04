#pragma once

class CNTparticle
{
private:
    int ID;
    
public:
    //Constructor and destructor
    CNTparticle();
    ~CNTparticle();
    
    void Initialize(int i);
    
    int GetID();
};
