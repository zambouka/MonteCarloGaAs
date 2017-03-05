#include <iostream>
#include "electron.h"

using namespace std;

Electron::Electron() //dummy code
{
    m = 9.11e-31;
    q = 1.60217e-19;
    h = 1.054571800e-34;
    valley = 0;
    effectiveMassG = 0.063;
    effectiveMassL = 0.17;
    effectiveMassX = 0.58;
}

Electron::~Electron()
{
    cout << "dead!" << endl;
    x = 0;
    y = 0;
    z = 0;
    energy = 0;
    kX = 0;
    kY = 0;
    kZ = 0;
    m = 0;
    q = 0;
    oldKX = 0;
    oldEnergy = 0;
}


void Electron::setKz(double _kZ)
{
    oldKZ = kZ;
    kZ = _kZ;
}
void Electron::setLocation(double _x, double _y, double _z)
{
    x = _x;
    y = _y;
    z = _z;
}

void Electron::setKx(double _kX)
{
    oldKX = kX;
    kX = _kX;
}

void Electron::setKy(double _kY)
{
    oldKY = kY;
    kY = _kY;
}

void Electron::setEnergy(double _energy)
{
    oldEnergy = energy;
    energy = _energy;
}

void Electron::setKDirection(double _kX, double _kY, double _kZ)
{
    kX = _kX;
    kY = _kY;
    kZ = _kZ;
}

double Electron::getOldKX()
{
    return oldKX;
}

double Electron::getOldKY()
{
    return oldKY;
}


double Electron::getkX()
{
    return kX;
}
double Electron::getkY()
{
    return kY;
}
double Electron::getkZ()
{
    return kZ;
}

double Electron::getX()
{
    return x;
}
double Electron::getY()
{
    return y;
}
double Electron::getZ()
{
    return z;
}
double Electron::getEnergy()
{
    return energy;
}

double Electron::getOldEnergy()
{
    return oldEnergy;
}



void Electron::recalculationOfEnergy()
{
    oldEnergy = energy;
    switch (valley)
    {
        case 0:
            energy = h*h*(kX*kX + kY*kY + kZ*kZ) / (2 * effectiveMassG*m *q);
            break;
        case 1:
            energy = h*h*(kX*kX + kY*kY + kZ*kZ) / (2 * effectiveMassL*m*q);
            break;
        case 2:
            energy = h*h*(kX*kX + kY*kY + kZ*kZ) / (2 * effectiveMassX*m*q);
            break;
    }
}

void Electron::setValley(int newValley)
{
    if ((0 <= newValley) && (newValley <= 2))
    {
        valley = newValley;
    }
    else
    {
        cout << "Wrong valley!";
    }
}

int Electron::getValley()
{
    return valley;
}

double Electron::getTau()
{
    return tau;
}

void Electron::setTau(double t)
{
    tau = t;
}





