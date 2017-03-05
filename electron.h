//
// Created by vanadiuz on 3/5/17.
//

#ifndef MCGAAS_ELECTRON_H
#define MCGAAS_ELECTRON_H

class Electron
{
public:
    Electron();
    ~Electron();
    void setKx(double _kX);
    void setKy(double _kY);
    void setKz(double _kZ);
    void setLocation(double _x, double _y, double _z);
    void setEnergy(double energy);
    void setKDirection(double _kX, double _kY, double _kZ);
    double getOldKX();
    double getOldKY();
    double getkX();
    double getkY();
    double getkZ();
    double getX();
    double getY();
    double getZ();
    double getEnergy();
    double getOldEnergy();
    void recalculationOfEnergy();
    void setValley(int);
    int getValley();
    double getTau();
    void setTau(double);

private:
    double x;
    double y;
    double z;
    double energy;
    double kX;
    double kY;
    double kZ;
    double m;
    double q;
    double oldKX;
    double oldKY;
    double oldKZ;
    double oldEnergy;
    int valley; //0 - G, 1 - L, 2 - X
    double effectiveMassG;
    double effectiveMassL;
    double effectiveMassX;
    double h;
    double tau;
};


#endif //MCGAAS_ELECTRON_H
