//
// Created by vanadiuz on 3/5/17.
//

#ifndef MCGAAS_MAINAREA_H
#define MCGAAS_MAINAREA_H

#include <vector>
class Electron;
class mainArea
{
public:
    mainArea(int);
    ~mainArea();



    void init(Electron &e);
    void process(Electron &e);
    void timeFreeMovement(Electron &e);
    void freeFlight(Electron &e);
    void drift(Electron &e, double t);
    void scatterCarrier(Electron &e);
    void analysis(Electron &e);
    void write();
    double driftVY(int n);
    void clearOut();

    //Gamma
    void makeTableTypesOfScatteringG();
    void selectTheTypeOfScatteringG(Electron &e);
    void changeValleyAngleLEnergyG(Electron &e, double);
    double changeEnergyG(int);
    double polarOpticalScatteringWithAbsorptionG(double energy);
    double polarOpticalScatteringWithEmissionG(double energy);
    double acousticScatteringG(double energy);
    double intervalleyScatteringNVGtoXAbsorbation(double energy);
    double intervalleyScatteringNVGtoXEmission(double energy);
    double intervalleyScatteringNVGtoLAbsorbation(double energy);
    double intervalleyScatteringNVGtoLEmission(double energy);
    double coulombG(double energy);
    //Gamma

    //L
    void makeTableTypesOfScatteringL();
    void selectTheTypeOfScatteringL(Electron &e);
    void changeValleyAngleLEnergyL(Electron &e, double);
    double changeEnergyL(int);
    double polarOpticalScatteringWithAbsorptionL(double energy);
    double polarOpticalScatteringWithEmissionL(double energy);
    double acousticScatteringL(double energy);
    double intervalleyScatteringNVLtoXAbsorbation(double energy);
    double intervalleyScatteringNVLtoXEmission(double energy);
    double intervalleyScatteringNVLtoGAbsorbation(double energy);
    double intervalleyScatteringNVLtoGEmission(double energy);
    double intervalleyScatteringEVLtoLEmission(double energy);
    double intervalleyScatteringEVLtoLAbsorbation(double energy);
    double coulombL(double energy);
    //L

    //X
    void makeTableTypesOfScatteringX();
    void selectTheTypeOfScatteringX(Electron &e);
    void changeValleyAngleLEnergyX(Electron &e, double);
    double changeEnergyX(int);
    double polarOpticalScatteringWithAbsorptionX(double energy);
    double polarOpticalScatteringWithEmissionX(double energy);
    double acousticScatteringX(double energy);
    double intervalleyScatteringNVXtoLAbsorbation(double energy);
    double intervalleyScatteringNVXtoLEmission(double energy);
    double intervalleyScatteringNVXtoGAbsorbation(double energy);
    double intervalleyScatteringNVXtoGEmission(double energy);
    double intervalleyScatteringEVXtoXEmission(double energy);
    double intervalleyScatteringEVXtoXAbsorbation(double energy);
    double coulombX(double energy);
    //X

    //angles
    void changeKDirectionPOSWE(Electron &e);
    void changeKDirectionPOSWA(Electron &e);
    void changeKDirectionAS(Electron &e);
    void changeKDirectionCoulombAngle(Electron &e);




private:
    double sigmaG;
    double sigmaL;
    double sigmaX;
    double pi;
    double phononGL;
    double phononGX;
    double phononLG;
    double phononLL;
    double phononLX;
    double phononXG;
    double phononXL;
    double phononXX;
    double eqValleysG;
    double eqValleysL;
    double eqValleysX;
    double q;
    double mq;
    double eps0;
    double movenentTime;
    double time;
    double newKX;
    double Ex;
    double h;
    double highFrecDielectrConst;
    double lowFrecDielectrConst;
    double LOPhononEnergy;
    double T;
    double kb;
    double acousticDefomationPotential;
    double density;
    double soundVelocity;
    double energySeparGL;
    double energySeparGX;
    double energySeparLX;
    double effectiveMassG;
    double effectiveMassL;
    double effectiveMassX;
    double deformationPotentialConstantGL;
    double deformationPotentialConstantGX;
    double deformationPotentialConstantLX;
    double deformationPotentialConstantLL;
    double deformationPotentialConstantXX;
    double intervalleyPhononEnergy;
    double nonparabolicityL;
    double nonparabolicityG;
    double nonparabolicityX;
    double equivalIntervalScatPhononEnergy;
    double doppingDensity;
    double tau;
    double totalTime;
    double fieldX;
    double fieldY;
    double fieldZ;

    double L;
    std::vector<std::vector<double>> *typesOfScatteringL;

    double X;
    std::vector<std::vector<double>> *typesOfScatteringX;

    double G;
    double G1;
    std::vector<std::vector<double>> *typesOfScatteringG;

    double velXSum[3];
    double velYSum[3];
    double velZSum[3];
    double energySum[3];
    double nValley[3];
    double velocityX[3];
    double velocityY[3];
    double velocityZ[3];
    //double waveVectorX[500];
    //double waveVectorY[500];
};




#endif //MCGAAS_MAINAREA_H
