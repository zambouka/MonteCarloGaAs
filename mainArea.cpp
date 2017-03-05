#include <iostream>
#include "mainArea.h"
#include <random>
#include <iomanip>
#include "electron.h"
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

mainArea::mainArea(int E)
{
    h = 1.054571800e-34;
    q = 1.6021766e-19;
    mq = 9.11e-31;
    eps0 = 8.85419e-12;
    highFrecDielectrConst = 10.92 * eps0;
    lowFrecDielectrConst = 12.9 * eps0;
    LOPhononEnergy = 0.03536;
    nonparabolicityL = 0.50;
    nonparabolicityG = 0.62;
    nonparabolicityX = 0.30;
    T = 400;
    kb = 1.38064852e-23;
    acousticDefomationPotential = 7;
    density = 5370;
    soundVelocity = 5.22e3;
    energySeparGL = 0.29;
    energySeparGX = 0.48;
    energySeparLX = energySeparGX- energySeparGL;
    eqValleysG = 1;
    eqValleysL = 4;
    eqValleysX = 3;
    effectiveMassG = 0.063;
    effectiveMassL = 0.17;// 0.17;
    effectiveMassX = 0.58;//0.58;
    deformationPotentialConstantGL = 0.18e11;
    deformationPotentialConstantGX = 1e11;
    deformationPotentialConstantLX = 0.1e11;
    deformationPotentialConstantLL = 0.5e11;
    deformationPotentialConstantXX = 1e11;
    phononGL = 0.0278;
    phononGX = 0.0299;
    phononLG = 0.0278;
    phononLL = 0.029;
    phononLX = 0.0293;
    phononXG = 0.0299;
    phononXL = 0.0293;
    phononXX = 0.0299;
    intervalleyPhononEnergy = 0.0299;
    equivalIntervalScatPhononEnergy = intervalleyPhononEnergy;
    pi = 3.14159265358979323846;
    sigmaG = 7.01;
    sigmaL = 9.2;
    sigmaX = 9.0;
    doppingDensity = 1e20;
    tau = 1e-14;
    totalTime = 100e-12;
    fieldX = 0;
    fieldY = E;
    fieldZ = 0;

    makeTableTypesOfScatteringG();
    makeTableTypesOfScatteringL();
    makeTableTypesOfScatteringX();

    for (size_t i = 0; i < 3; i++)
    {
        velXSum[i] = 0;
        velYSum[i] = 0;
        velZSum[i] = 0;
        energySum[i] = 0;
        nValley[i] = 0;
    }
    //for (size_t i = 0; i < 500; i++)
    //{
    //	waveVectorX[i] = 0;
    //	waveVectorY[i] = 0;
    //}
}


mainArea::~mainArea()
{
}

void mainArea::init(Electron & e)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distr(0, 1);
    double r = distr(gen);
    double energy = -(1.5*kb*T / q)*log(r);
    e.setEnergy(energy);
    //energySum[0] += energy;

    e.setValley(0);

    double k = sqrt(2 * mq*effectiveMassG)*sqrt(q) / h * sqrt(e.getEnergy()*(1+nonparabolicityG*e.getEnergy()));
    double fai = 2 * pi*r;
    double ct = 1 - 2 * r;
    double st = sqrt(1 - ct*ct);

    e.setKx(k*st*cos(fai));
    e.setKy(k*st*sin(fai));
    e.setKz(k*ct);

    timeFreeMovement(e);
    e.setTau(movenentTime);
}

void mainArea::process(Electron & e)
{
    freeFlight(e);
}

void mainArea::timeFreeMovement(Electron &e)
{
    random_device rd;
    mt19937 gen(rd());
    double max = 0;
    switch (e.getValley())
    {
        case 0:
            if (e.getEnergy() <= 0.5)
            {
                max = G1;
            }
            else
            {
                max = G;
            }
            break;
        case 1:
            max = L;
            break;
        case 2:
            max = X;
            break;
    }
    //exponential_distribution<> distr(max);
    uniform_real_distribution<> d(0, 1);
    double r = d(gen);
    double gt = -log(r)*1/max;
    while (r <= 1.e-6)
    {
        r = d(gen);
        gt = -log(r) * 1 / max;
    }
    movenentTime = gt;
}

void mainArea::freeFlight(Electron & e)
{
    double dte = e.getTau();
    double dt2;
    if (dte >= tau)
    {
        dt2 = tau;
    }
    else
    {
        dt2 = dte;
    }
    drift(e, dt2);
    while (dte < tau)
    {
        double dte2 = dte;
        scatterCarrier(e);
        timeFreeMovement(e);
        double dt3 = movenentTime;
        double dtp = tau - dte2;
        if (dt3 <= dtp)
        {
            dt2 = dt3;
        }
        else
        {
            dt2 = dtp;
        }
        drift(e, dt2);
        dte2 = dte2 + dt3;
        dte = dte2;
    }
    if (dte > tau)
    {
        e.setTau(dte - tau);
    }
}

void mainArea::drift(Electron & e, double t)
{
    double qhl = q/h*t;
    double dkx = -qhl*fieldX;
    double dky = -qhl*fieldY;
    double dkz = -qhl*fieldZ;

    e.setKx(e.getkX() + dkx);
    e.setKy(e.getkY() + dky);
    e.setKz(e.getkZ() + dkz);

    double skx = e.getkX()*e.getkX();
    double sky = e.getkY()*e.getkY();
    double skz = e.getkZ()*e.getkZ();
    double sk = skx + sky + skz;
    double k = sqrt(sk);
    double m;
    double alpha;
    switch (e.getValley())
    {
        case 0:
            m = effectiveMassG *mq;
            alpha = nonparabolicityG;
            break;
        case 1:
            m = effectiveMassL *mq;
            alpha = nonparabolicityL;
            break;
        case 2:
            m = effectiveMassX *mq;
            alpha = nonparabolicityX;
            break;
    }
    double gk = h*h / (2 * m*q)*sk;
    double energy = 2 * gk / (1 + sqrt(1 + 4 * alpha*gk));
    e.setEnergy(energy);
}

void mainArea::scatterCarrier(Electron & e)
{
    if (e.getValley() == 0)
    {
        selectTheTypeOfScatteringG(e);
    }
    if (e.getValley() == 1)
    {
        selectTheTypeOfScatteringL(e);
    }
    if (e.getValley() == 2)
    {
        selectTheTypeOfScatteringX(e);
    }
}

void mainArea::analysis(Electron & e)
{

    double v = e.getValley();
    double energy = e.getEnergy();
    double alpha;
    double m;
    switch (e.getValley())
    {
        case 0:
            m = mq*effectiveMassG;
            alpha = nonparabolicityG;
            break;
        case 1:
            m = mq*effectiveMassL;
            alpha = nonparabolicityL;
            break;
        case 2:
            m = mq*effectiveMassX;
            alpha = nonparabolicityX;
            break;
    }
    double denom = 1 / (1 + 2 * energy*alpha);
    double velX = h*e.getkX()*denom / m;
    double velY = h*e.getkY()*denom / m;
    double velZ = h*e.getkZ()*denom / m;
    velXSum[e.getValley()] += velX;
    velYSum[e.getValley()] += velY;
    velZSum[e.getValley()] += velZ;
    energySum[e.getValley()] += energy;
    nValley[e.getValley()] += 1;

    //if (abs(e.getkX()) <= 1.5e9)
    //{
    //	double wv = e.getkX() + 1.5e9;
    //	int n = (wv / 3e9)*500;
    //	waveVectorX[n] ++;
    //}
    //if (abs(e.getkY()) <= 1.5e9)
    //{
    //	double wv = e.getkY() + 1.5e9;
    //	int n = (wv / 3e9) * 500;
    //	waveVectorY[n] ++;
    //}
}

void mainArea::write()
{
    ofstream vel("velocity.txt");
    ofstream time("time.txt");
    ofstream energy("energy.txt");
    ofstream oWaveVectorX("kX.txt");
    ofstream oWaveVectorY("kY.txt");
    for (size_t i = 0; i < 3; i++)
    {
        if (nValley[i] != 0)
        {
            velocityX[i] = velXSum[i] / nValley[i];
            velocityY[i] = velYSum[i] / nValley[i];
            velocityZ[i] = velZSum[i] / nValley[i];
            energySum[i] = energySum[i] / nValley[i];
            vel << velocityX[i] << "      " << velocityY[i] << "      " << velocityZ[i] << endl;
            energy << energySum[i] << endl;
        }
    }
    //for (size_t i = 0; i < 500; i++)
    //{
    //	oWaveVectorX << waveVectorX[i] << endl;
    //	oWaveVectorY << waveVectorY[i] << endl;
    //}
}

double mainArea::driftVY(int n)
{
    double v = 0;
    for (size_t i = 0; i < 3; i++)
    {
        v += velocityY[i] * nValley[i] / n;
    }
    return v;
}

void mainArea::clearOut()
{
    for (size_t i = 0; i < 3; i++)
    {
        velXSum[i] = 0;
        velYSum[i] = 0;
        velZSum[i] = 0;
        energySum[i] = 0;
        nValley[i] = 0;
    }
}



void mainArea::makeTableTypesOfScatteringG() // 0 - 0,5 eV .... step = (0.005)
{
    typesOfScatteringG = new vector<vector<double>>(600, vector<double>(9));
    G = 0;
    G1 = 0;
    for (size_t i = 1; i < 600; i++)
    {
        (*typesOfScatteringG)[i][0] = polarOpticalScatteringWithAbsorptionG(0.005*(i));
        (*typesOfScatteringG)[i][1] = polarOpticalScatteringWithEmissionG(0.005*(i));
        (*typesOfScatteringG)[i][2] = acousticScatteringG(0.005*(i));
        (*typesOfScatteringG)[i][3] = intervalleyScatteringNVGtoXAbsorbation(0.005*(i));
        (*typesOfScatteringG)[i][4] = intervalleyScatteringNVGtoXEmission(0.005*(i));
        (*typesOfScatteringG)[i][5] = intervalleyScatteringNVGtoLAbsorbation(0.005*(i));
        (*typesOfScatteringG)[i][6] = intervalleyScatteringNVGtoLEmission(0.005*(i));
        (*typesOfScatteringG)[i][7] = coulombG(0.005*(i));

        double gg = 0;
        for (size_t j = 0; j < 8; j++)
        {
            gg += (*typesOfScatteringG)[i][j];
        }
        if (gg > G)
        {
            G = gg;
        }
        if (i < 100)
        {
            if (gg > G1)
            {
                G1 = gg;
            }
        }
    }
    for (size_t i = 100; i < 600; i++)
    {
        (*typesOfScatteringG)[i][8] = G;
    }
    for (size_t i = 1; i < 100; i++)
    {
        (*typesOfScatteringG)[i][8] = G1;
    }
}

void mainArea::selectTheTypeOfScatteringG(Electron &e)
{
    random_device rd;
    mt19937 gen(rd());
    double max;
    if (e.getEnergy() <= 0.5)
    {
        max = G1;
    }
    else
    {
        max = G;
    }
    uniform_real_distribution<> distr(0, max);
    double r = distr(gen);
    double sum = 0;
    int typeOfScattering;
    int type = e.getEnergy() / 0.005;
    if (type > 599)
    {
        type = 599;
    }
    if (e.getEnergy() < 0.005)
    {
        type = 1;
    }

    for (size_t i = 0; i < 9; i++)
    {
        sum += (*typesOfScatteringG)[type][i];
        if (r < sum)
        {
            typeOfScattering = i;
            break;
        }
    }



    if (typeOfScattering != 8)
    {
        changeValleyAngleLEnergyG(e, typeOfScattering);
    }
}

void mainArea::changeValleyAngleLEnergyG(Electron & e, double typeOfScattering)
{
    if ((typeOfScattering == 3) || (typeOfScattering == 4))
    {
        e.setValley(2);
    }
    if ((typeOfScattering == 5) || (typeOfScattering == 6))
    {
        e.setValley(1);
    }

    e.setEnergy(e.getEnergy() + changeEnergyG(typeOfScattering));

    if (typeOfScattering == 0)
    {
        changeKDirectionPOSWA(e);
    }
    if (typeOfScattering == 1)
    {
        changeKDirectionPOSWE(e);
    }
    if ((typeOfScattering >= 2) && (typeOfScattering <= 6))
    {
        changeKDirectionAS(e);
    }
    if (typeOfScattering == 7)
    {
        changeKDirectionCoulombAngle(e);
    }
}

double mainArea::changeEnergyG(int type)
{
    if (type == 0)
    {
        return LOPhononEnergy;
    }
    if (type == 1)
    {
        return -LOPhononEnergy;
    }
    if (type == 3)
    {
        return phononGX - energySeparGX;
    }
    if (type == 4)
    {
        return -phononGX - energySeparGX;
    }
    if (type == 5)
    {
        return  phononGL - energySeparGL;
    }
    if (type == 6)
    {
        return -phononGL - energySeparGL;
    }
    return 0;
}

double mainArea::polarOpticalScatteringWithAbsorptionG(double energy)
{
    double alpha = nonparabolicityG;
    double N = (1 / (exp(LOPhononEnergy / (kb*T / q)) - 1));
    double c = pow(q, 3)*LOPhononEnergy*sqrt((effectiveMassG*mq) / (2 * q)) / (4 * pi * pow(h, 2))*((1 / highFrecDielectrConst) - (1 / lowFrecDielectrConst));
    double newEnergy = energy + LOPhononEnergy;
    double lamb = energy*(1 + alpha*energy);
    double lambNew = newEnergy*(1 + alpha*newEnergy);
    double A = pow((2 * (1 + 2 * alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew)), 2);
    double B = -2 * alpha*sqrt(lamb)*sqrt(lambNew)*(4 * (1 + alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew));
    double C = 4 * (1 + alpha*energy)*(1 + alpha*newEnergy)*(1 + 2 * alpha*energy)*(1 + 2 * alpha*newEnergy);
    double ab = abs((sqrt(lamb) + sqrt(lambNew)) / (sqrt(lamb) - sqrt(lambNew)));
    A = C = 4; B = 0; //!!!!
    double factor = (1 + 2 * alpha*newEnergy) / sqrt(lamb)*(A*log(ab) + B) / C;
    return N*c*factor;
}

double mainArea::polarOpticalScatteringWithEmissionG(double energy)
{
    if (energy < LOPhononEnergy)
    {
        return 0;
    }
    double alpha = nonparabolicityG;
    double N = (1 / (exp(LOPhononEnergy / (kb*T / q)) - 1) + 1);
    double c = pow(q, 3)*LOPhononEnergy*sqrt(effectiveMassG*mq / (2 * q)) / (4 * pi*pow(h, 2))*((1 / highFrecDielectrConst) - (1 / lowFrecDielectrConst));
    double newEnergy = energy - LOPhononEnergy;
    double lamb = energy*(1 + alpha*energy);
    double lambNew = newEnergy*(1 + alpha*newEnergy);
    double A = pow((2 * (1 + 2 * alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew)), 2);
    double B = -2 * alpha*sqrt(lamb)*sqrt(lambNew)*(4 * (1 + alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew));
    double C = 4 * (1 + alpha*energy)*(1 + alpha*newEnergy)*(1 + 2 * alpha*energy)*(1 + 2 * alpha*newEnergy);
    double ab = abs((sqrt(lamb) + sqrt(lambNew)) / (sqrt(lamb) - sqrt(lambNew)));
    A = C = 4; B = 0; //!!!!
    double factor = (1 + 2 * alpha*newEnergy) / sqrt(lamb)*(A*log(ab) + B) / C;
    double a = N*c*factor;
    return a;
}



void mainArea::changeKDirectionPOSWE(Electron &e)
{
    double m;
    double alpha;
    switch (e.getValley())
    {
        case 0:
            m = mq*effectiveMassG;
            alpha = nonparabolicityG;
            break;
        case 1:
            m = mq*effectiveMassL;
            alpha = nonparabolicityL;
            break;
        case 2:
            m = mq*effectiveMassX;
            alpha = nonparabolicityX;
            break;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> d(0, 1);

    double enew = e.getEnergy();
    double eold = e.getEnergy() + LOPhononEnergy;
    double kxy = sqrt(e.getkX()*e.getkX() + e.getkY()*e.getkY());
    double k = sqrt(kxy*kxy + e.getkZ()*e.getkZ());
    double cth0 = e.getkZ() / k;
    double sth0 = kxy / k;
    double cfi0 = e.getkX() / kxy;
    double sfi0 = e.getkY() / kxy;

    double kp = sqrt(2 * m*q) / h*sqrt(enew * (1 + alpha*enew));
    double ge = eold * (1 + alpha * eold);
    double gnew = enew * (1 + alpha*enew);
    double zeta = 2 * sqrt(ge*gnew) / (ge + gnew - 2 * sqrt(ge*gnew));
    double rr = d(gen);
    double cth = ((zeta + 1) - pow((2 * zeta + 1), rr)) / zeta;
    double sth = sqrt(1 - cth*cth);
    double fi = 2 * pi*d(gen);
    double cfi = cos(fi);
    double sfi = sin(fi);
    double kxp = kp*sth*cfi;
    double kyp = kp*sth*sfi;
    double kzp = kp*cth;

    double kx = kxp*cfi0*cth0 - kyp*sfi0 + kzp*cfi0*sth0;
    double ky = kxp*sfi0*cth0 + kyp*cfi0 + kzp*sfi0*sth0;
    double kz = -kxp*sth0 + kzp*cth0;

    e.setKDirection(kx, ky, kz);
}

void mainArea::changeKDirectionPOSWA(Electron &e)
{
    double m;
    double alpha;
    switch (e.getValley())
    {
        case 0:
            m = mq*effectiveMassG;
            alpha = nonparabolicityG;
            break;
        case 1:
            m = mq*effectiveMassL;
            alpha = nonparabolicityL;
            break;
        case 2:
            m = mq*effectiveMassX;
            alpha = nonparabolicityX;
            break;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> d(0, 1);

    double enew = e.getEnergy();
    double eold = e.getEnergy() - LOPhononEnergy;
    double kxy = sqrt(e.getkX()*e.getkX() + e.getkY()*e.getkY());
    double k = sqrt(kxy*kxy + e.getkZ()*e.getkZ());
    double cth0 = e.getkZ() / k;
    double sth0 = kxy / k;
    double cfi0 = e.getkX() / kxy;
    double sfi0 = e.getkY() / kxy;

    double kp = sqrt(2 * m*q) / h*sqrt(enew * (1 + alpha*enew));
    double ge = eold * (1 + alpha * eold);
    double gnew = enew * (1 + alpha*enew);
    double zeta = 2 * sqrt(ge*gnew) / (ge + gnew - 2 * sqrt(ge*gnew));
    double rr = d(gen);
    double cth = ((zeta + 1) - pow((2 * zeta + 1), rr)) / zeta;
    double sth = sqrt(1 - cth*cth);
    double fi = 2 * pi*d(gen);
    double cfi = cos(fi);
    double sfi = sin(fi);
    double kxp = kp*sth*cfi;
    double kyp = kp*sth*sfi;
    double kzp = kp*cth;

    double kx = kxp*cfi0*cth0 - kyp*sfi0 + kzp*cfi0*sth0;
    double ky = kxp*sfi0*cth0 + kyp*cfi0 + kzp*sfi0*sth0;
    double kz = -kxp*sth0 + kzp*cth0;

    e.setKDirection(kx, ky, kz);
}

double mainArea::acousticScatteringG(double energy)
{
    double lamb = energy*(1 + nonparabolicityG*energy);
    double j = (sqrt(2)*pow((effectiveMassG*mq), 1.5) * pow(q, 2.5) * kb * T * pow(sigmaG, 2)) / (pow(soundVelocity, 2)*density*pi*pow(h, 4))*sqrt(lamb)*(1 + 2 * nonparabolicityG*energy);
    return j;
}

void mainArea::changeKDirectionAS(Electron & e)
{
    double m;
    double alpha;
    switch (e.getValley())
    {
        case 0:
            m = mq*effectiveMassG;
            alpha = nonparabolicityG;
            break;
        case 1:
            m = mq*effectiveMassL;
            alpha = nonparabolicityL;
            break;
        case 2:
            m = mq*effectiveMassX;
            alpha = nonparabolicityX;
            break;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> d(0, 1);

    double rknew =sqrt(2*m*q)*sqrt(e.getEnergy()*(1 + alpha*e.getEnergy()))/h;
    double fi = 2 * pi *d(gen);
    double ct = 1-2 * d(gen);
    double st = sqrt(1 - ct*ct);
    double kx = rknew*st*cos(fi);
    double ky = rknew*st*sin(fi);
    double kz = rknew*ct;

    e.setKDirection(kx, ky, kz);
}

double mainArea::intervalleyScatteringNVGtoXAbsorbation(double energy)
{
    if (energy < (energySeparGX - phononGX))
    {
        return 0.0;
    }
    double N = 1 / (exp(phononGX / (kb*T / q)) - 1);
    double finalMass = effectiveMassX*mq;
    double c = eqValleysX*pow(deformationPotentialConstantGX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononGX)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononGX - energySeparGX;
    double lamb = ef*(1 + nonparabolicityX*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityX*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVGtoXEmission(double energy)
{
    if (energy < (energySeparGX + phononGX))
    {
        return 0.0;
    }
    double N = 1 / (exp(phononGX / (kb*T / q)) - 1) + 1;
    double finalMass = effectiveMassX*mq;
    double c = eqValleysX*pow(deformationPotentialConstantGX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononGX)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononGX - energySeparGX;
    double lamb = ef*(1 + nonparabolicityX*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityX*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVGtoLAbsorbation(double energy)
{
    if (energy < (energySeparGL - phononGL))
    {
        return 0.0;
    }
    double N = 1 / (exp(phononGL / (kb*T / q)) - 1);
    double finalMass = effectiveMassL*mq;
    double c = eqValleysL*pow(deformationPotentialConstantGL, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononGL)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononGL - energySeparGL;
    double lamb = ef*(1 + nonparabolicityL*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityL*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVGtoLEmission(double energy)
{
    if (energy < (energySeparGL + phononGL))
    {
        return 0.0;
    }
    double N = 1 / (exp(phononGL / (kb*T / q)) - 1) + 1;
    double finalMass = effectiveMassL*mq;
    double c = eqValleysL*pow(deformationPotentialConstantGL, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononGL)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononGL - energySeparGL;
    double lamb = ef*(1 + nonparabolicityL*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityL*ef);
    return ab*factor;
}

double mainArea::coulombG(double energy)
{
    double alpha = nonparabolicityG;
    double m = mq*effectiveMassG;
    double ld = sqrt(highFrecDielectrConst * kb*T / q * 1 / q * 1 / doppingDensity);
    double ed = h*h / (m * q * 2 * ld*ld);
    double c = doppingDensity*m*q / h*sqrt(2 * m) * 1 / pi * q / h * pow((ld*ld*q / (highFrecDielectrConst*h)), 2) *  sqrt(q);
    double a = energy*(1 + alpha*energy);
    double factor = sqrt(a) * (1 + 2 * alpha*energy) * 1 / (1 + 4 * a / ed);
    return c*factor;
}

void mainArea::makeTableTypesOfScatteringL()
{
    typesOfScatteringL = new vector<vector<double>>(600, vector<double>(11));
    L = 0;
    for (size_t i = 1; i < 600; i++)
    {
        (*typesOfScatteringL)[i][0] = polarOpticalScatteringWithAbsorptionL(0.005*(i));
        (*typesOfScatteringL)[i][1] = polarOpticalScatteringWithEmissionL(0.005*(i));
        (*typesOfScatteringL)[i][2] = acousticScatteringL(0.005*(i));
        (*typesOfScatteringL)[i][3] = intervalleyScatteringNVLtoXAbsorbation(0.005*(i));
        (*typesOfScatteringL)[i][4] = intervalleyScatteringNVLtoXEmission(0.005*(i));
        (*typesOfScatteringL)[i][5] = intervalleyScatteringNVLtoGAbsorbation(0.005*(i));
        (*typesOfScatteringL)[i][6] = intervalleyScatteringNVLtoGEmission(0.005*(i));
        (*typesOfScatteringL)[i][7] = intervalleyScatteringEVLtoLAbsorbation(0.005*(i));
        (*typesOfScatteringL)[i][8] = intervalleyScatteringEVLtoLEmission(0.005*(i));
        (*typesOfScatteringL)[i][9] = coulombL(0.005*(i));

        double ll = 0;
        for (size_t j = 0; j < 10; j++)
        {
            ll += (*typesOfScatteringL)[i][j];
        }
        if (ll > L)
        {
            L = ll;
        }
    }
    for (size_t i = 1; i < 600; i++)
    {
        (*typesOfScatteringL)[i][10] = L;
    }
}

void mainArea::selectTheTypeOfScatteringL(Electron & e)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distr(0, L);
    double r = distr(gen);
    double sum = 0;
    int type = e.getEnergy() / 0.005;
    int typeOfScattering;
    if (type > 599)
    {
        type = 599;
    }
    if (e.getEnergy() < 0.005)
    {
        type = 1;
    }

    for (size_t i = 0; i < 11; i++)
    {
        sum += (*typesOfScatteringL)[type][i];
        if (r <= sum)
        {
            typeOfScattering = i;
            break;
        }
    }

    if (typeOfScattering != 10)
    {
        changeValleyAngleLEnergyL(e, typeOfScattering);
    }
}

void mainArea::changeValleyAngleLEnergyL(Electron & e, double typeOfScattering)
{
    if ((typeOfScattering == 3) || (typeOfScattering == 4))
    {
        e.setValley(2);
    }
    if ((typeOfScattering == 5) || (typeOfScattering == 6))
    {
        e.setValley(0);
    }

    e.setEnergy(e.getEnergy() + changeEnergyL(typeOfScattering));

    if (typeOfScattering == 0)
    {
        changeKDirectionPOSWA(e);
    }
    if (typeOfScattering == 1)
    {
        changeKDirectionPOSWE(e);
    }
    if ((typeOfScattering >= 2) && (typeOfScattering <= 8))
    {
        changeKDirectionAS(e);
    }
    if (typeOfScattering == 9)
    {
        changeKDirectionCoulombAngle(e);
    }
}

double mainArea::changeEnergyL(int type)
{
    if (type == 0)
    {
        return LOPhononEnergy;
    }
    if (type == 1)
    {
        return -LOPhononEnergy;
    }
    if (type == 3)
    {
        return -(energySeparLX - phononLX);
    }
    if (type == 4)
    {
        return -(energySeparLX + phononLX);
    }
    if (type == 5)
    {
        return  phononLG + energySeparGL;
    }
    if (type == 6)
    {
        return  -phononLG + energySeparGL;
    }
    if (type == 7)
    {
        return  phononLL;
    }
    if (type == 8)
    {
        return  -phononLL;
    }
    return 0;
}

double mainArea::polarOpticalScatteringWithAbsorptionL(double energy)
{
    double alpha = nonparabolicityL;
    double N = (1 / (exp(LOPhononEnergy / (kb*T / q)) - 1));
    double c = pow(q, 3)*LOPhononEnergy*sqrt((effectiveMassL*mq) / (2 * q)) / (4 * pi * pow(h,2))*((1 / highFrecDielectrConst) - (1 / lowFrecDielectrConst));
    double newEnergy = energy + LOPhononEnergy;
    double lamb = energy*(1 + alpha*energy);
    double lambNew = newEnergy*(1 + alpha*newEnergy);
    double A = pow((2 * (1 + 2 * alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew)), 2);
    double B = -2 * alpha*sqrt(lamb)*sqrt(lambNew)*(4 * (1 + alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew));
    double C = 4 * (1 + alpha*energy)*(1 + alpha*newEnergy)*(1 + 2 * alpha*energy)*(1 + 2 * alpha*newEnergy);
    double ab = abs((sqrt(lamb) + sqrt(lambNew)) / (sqrt(lamb) - sqrt(lambNew)));
    A = C = 4; B = 0; //!!!!
    double factor = (1 + 2 * alpha*newEnergy) / sqrt(lamb)*(A*log(ab) + B) / C;
    return N*c*factor;
}

double mainArea::polarOpticalScatteringWithEmissionL(double energy)
{
    if (energy < LOPhononEnergy)
    {
        return 0;
    }
    double alpha = nonparabolicityL;
    double N = (1 / (exp(LOPhononEnergy / (kb*T/q)) - 1) + 1);
    double c = pow(q, 3)*LOPhononEnergy*sqrt(effectiveMassL*mq / (2 * q)) / (4 * pi*pow(h, 2))*((1 / highFrecDielectrConst) - (1 / lowFrecDielectrConst));
    double newEnergy = energy - LOPhononEnergy;
    double lamb = energy*(1 + alpha*energy);
    double lambNew = newEnergy*(1 + alpha*newEnergy);
    double A = pow((2 * (1 + 2 * alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew)), 2);
    double B = -2 * alpha*sqrt(lamb)*sqrt(lambNew)*(4 * (1 + alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew));
    double C = 4 * (1 + alpha*energy)*(1 + alpha*newEnergy)*(1 + 2 * alpha*energy)*(1 + 2 * alpha*newEnergy);
    double ab = abs((sqrt(lamb) + sqrt(lambNew)) / (sqrt(lamb) - sqrt(lambNew)));
    A = C = 4; B = 0; //!!!!
    double factor = (1 + 2 * alpha*newEnergy) / sqrt(lamb)*(A*log(ab) + B) / C;
    double a = N*c*factor;
    return a;
}

double mainArea::acousticScatteringL(double energy)
{
    double lamb = energy*(1 + nonparabolicityL*energy);
    double j = (sqrt(2)*pow((effectiveMassL*mq), 1.5) * pow(q, 2.5) * kb * T * pow(sigmaL, 2)) / (pow(soundVelocity, 2)*density*pi*pow(h, 4))*sqrt(lamb)*(1 + 2 * nonparabolicityL*energy);
    return j;
}

double mainArea::intervalleyScatteringNVLtoXAbsorbation(double energy)
{
    if (energy < (energySeparLX - phononLX))
    {
        return 0.0;
    }
    double N = 1 / (exp(phononLX / (kb*T / q)) - 1);
    double finalMass = effectiveMassX*mq;
    double c = eqValleysX*pow(deformationPotentialConstantLX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononLX)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononLX - energySeparLX;
    double lamb = ef*(1 + nonparabolicityX*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityX*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVLtoXEmission(double energy)
{
    if (energy < (energySeparLX + phononLX))
    {
        return 0.0;
    }
    double N = 1 / (exp(phononLX / (kb*T / q)) - 1) + 1;
    double finalMass = effectiveMassX*mq;
    double c = eqValleysX*pow(deformationPotentialConstantLX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononLX)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononLX - energySeparLX;
    double lamb = ef*(1 + nonparabolicityX*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityX*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVLtoGAbsorbation(double energy)
{
    double N = 1 / (exp(phononLG / (kb*T/q)) - 1);
    double finalMass = effectiveMassG*mq;
    double c = eqValleysG*pow(deformationPotentialConstantGL, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononLG)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononLG + energySeparGL;
    double lamb = ef*(1 + nonparabolicityG*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityG*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVLtoGEmission(double energy)
{
    double N = 1 / (exp(phononLG / (kb*T / q)) - 1) + 1;
    double finalMass = effectiveMassG*mq;
    double c = eqValleysG*pow(deformationPotentialConstantGL, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononLG)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononLG + energySeparGL;
    double lamb = ef*(1 + nonparabolicityG*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityG*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringEVLtoLEmission(double energy)
{
    if (energy < phononLL)
    {
        return 0.0;
    }
    double N = 1 / (exp(phononLL / (kb*T / q)) - 1) + 1;
    double finalMass = effectiveMassL*mq;
    double c = (eqValleysL-1)*pow(deformationPotentialConstantLL, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononLL)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononLL;
    double lamb = ef*(1 + nonparabolicityL*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityL*ef);
    return ab*factor;

}

double mainArea::intervalleyScatteringEVLtoLAbsorbation(double energy)
{
    double N = 1 / (exp(phononLL / (kb*T / q)) - 1);
    double finalMass = effectiveMassL*mq;
    double c = (eqValleysL - 1)*pow(deformationPotentialConstantLL, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononLL)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononLL;
    double lamb = ef*(1 + nonparabolicityL*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityL*ef);
    return ab*factor;

}

double mainArea::coulombL(double energy)
{
    double alpha = nonparabolicityL;
    double m = mq*effectiveMassL;
    double ld = sqrt(highFrecDielectrConst * kb*T / q * 1 / q * 1 / doppingDensity);
    double ed = h*h / (m * q * 2 * ld*ld);
    double c = doppingDensity*m*q / h*sqrt(2 * m) * 1 / pi * q / h * pow((ld*ld*q / (highFrecDielectrConst*h)), 2) * sqrt(q);
    double a = energy*(1 + alpha*energy);
    double factor = sqrt(a) * (1 + 2 * alpha*energy) * 1 / (1 + 4 * a / ed);
    return c*factor;
}

void mainArea::makeTableTypesOfScatteringX()
{
    typesOfScatteringX = new vector<vector<double>>(600, vector<double>(11));
    X = 0;
    for (size_t i = 1; i < 600; i++)
    {
        (*typesOfScatteringX)[i][0] = polarOpticalScatteringWithAbsorptionX(0.005*(i));
        (*typesOfScatteringX)[i][1] = polarOpticalScatteringWithEmissionX(0.005*(i));
        (*typesOfScatteringX)[i][2] = acousticScatteringX(0.005*(i));
        (*typesOfScatteringX)[i][3] = intervalleyScatteringNVXtoLAbsorbation(0.005*(i));
        (*typesOfScatteringX)[i][4] = intervalleyScatteringNVXtoLEmission(0.005*(i));
        (*typesOfScatteringX)[i][5] = intervalleyScatteringNVXtoGAbsorbation(0.005*(i));
        (*typesOfScatteringX)[i][6] = intervalleyScatteringNVXtoGEmission(0.005*(i));
        (*typesOfScatteringX)[i][7] = intervalleyScatteringEVXtoXAbsorbation(0.005*(i));
        (*typesOfScatteringX)[i][8] = intervalleyScatteringEVXtoXEmission(0.005*(i));
        (*typesOfScatteringX)[i][9] = coulombX(0.005*(i));

        double xx = 0;
        for (size_t j = 0; j < 10; j++)
        {
            xx += (*typesOfScatteringX)[i][j];
        }
        if (xx > X)
        {
            X = xx;
        }
    }
    for (size_t i = 1; i < 600; i++)
    {
        (*typesOfScatteringX)[i][10] = X;
    }
}

void mainArea::selectTheTypeOfScatteringX(Electron & e)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distr(0, X);
    double r = distr(gen);
    double sum = 0;
    int type = e.getEnergy() / 0.005;
    int typeOfScattering;
    if (type > 599)
    {
        type = 599;
    }
    if (e.getEnergy() < 0.005)
    {
        type = 1;
    }

    for (size_t i = 0; i < 11; i++)
    {
        sum += (*typesOfScatteringX)[type][i];
        if (r <= sum)
        {
            typeOfScattering = i;
            break;
        }
    }
    if (typeOfScattering != 10)
    {
        changeValleyAngleLEnergyX(e, typeOfScattering);
    }
}

void mainArea::changeValleyAngleLEnergyX(Electron & e, double typeOfScattering)
{
    if ((typeOfScattering == 3) || (typeOfScattering == 4))
    {
        e.setValley(1);
    }
    if ((typeOfScattering == 5) || (typeOfScattering == 6))
    {
        e.setValley(0);
    }

    e.setEnergy(e.getEnergy() + changeEnergyX(typeOfScattering));

    if (typeOfScattering == 0)
    {
        changeKDirectionPOSWA(e);
    }
    if (typeOfScattering == 1)
    {
        changeKDirectionPOSWE(e);
    }
    if ((typeOfScattering >= 2) && (typeOfScattering <= 8))
    {
        changeKDirectionAS(e);
    }
    if (typeOfScattering == 9)
    {
        changeKDirectionCoulombAngle(e);
    }
}

double mainArea::changeEnergyX(int type)
{
    if (type == 0)
    {
        return LOPhononEnergy;
    }
    if (type == 1)
    {
        return -LOPhononEnergy;
    }
    if (type == 3)
    {
        return energySeparLX + phononLX;
    }
    if (type == 4)
    {
        return energySeparLX - phononLX;
    }
    if (type == 5)
    {
        return  phononXG + energySeparGX;
    }
    if (type == 6)
    {
        return  -phononXG + energySeparGX;
    }
    if (type == 7)
    {
        return  phononXX;
    }
    if (type == 8)
    {
        return  -phononXX;
    }
    return 0;
}

double mainArea::polarOpticalScatteringWithAbsorptionX(double energy)
{
    double alpha = nonparabolicityX;
    double N = (1 / (exp(LOPhononEnergy / (kb*T / q)) - 1));
    double c = pow(q, 3)*LOPhononEnergy*sqrt((effectiveMassX*mq) / (2 * q)) / (4 * pi * pow(h, 2))*((1 / highFrecDielectrConst) - (1 / lowFrecDielectrConst));
    double newEnergy = energy + LOPhononEnergy;
    double lamb = energy*(1 + alpha*energy);
    double lambNew = newEnergy*(1 + alpha*newEnergy);
    double A = pow((2 * (1 + 2 * alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew)), 2);
    double B = -2 * alpha*sqrt(lamb)*sqrt(lambNew)*(4 * (1 + alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew));
    double C = 4 * (1 + alpha*energy)*(1 + alpha*newEnergy)*(1 + 2 * alpha*energy)*(1 + 2 * alpha*newEnergy);
    double ab = abs((sqrt(lamb) + sqrt(lambNew)) / (sqrt(lamb) - sqrt(lambNew)));
    A = C = 4; B = 0; //!!!!
    double factor = (1 + 2 * alpha*newEnergy) / sqrt(lamb)*(A*log(ab) + B) / C;
    return N*c*factor;
}

double mainArea::polarOpticalScatteringWithEmissionX(double energy)
{
    if (energy < LOPhononEnergy)
    {
        return 0;
    }
    double alpha = nonparabolicityX;
    double N = (1 / (exp(LOPhononEnergy / (kb*T / q)) - 1))+1;
    double c = pow(q, 3)*LOPhononEnergy*sqrt((effectiveMassX*mq) / (2 * q)) / (4 * pi * pow(h, 2))*((1 / highFrecDielectrConst) - (1 / lowFrecDielectrConst));
    double newEnergy = energy - LOPhononEnergy;
    double lamb = energy*(1 + alpha*energy);
    double lambNew = newEnergy*(1 + alpha*newEnergy);
    double A = pow((2 * (1 + 2 * alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew)), 2);
    double B = -2 * alpha*sqrt(lamb)*sqrt(lambNew)*(4 * (1 + alpha*energy)*(1 + alpha*newEnergy) + alpha*(lamb + lambNew));
    double C = 4 * (1 + alpha*energy)*(1 + alpha*newEnergy)*(1 + 2 * alpha*energy)*(1 + 2 * alpha*newEnergy);
    double ab = abs((sqrt(lamb) + sqrt(lambNew)) / (sqrt(lamb) - sqrt(lambNew)));
    A = C = 4; B = 0; //!!!!
    double factor = (1 + 2 * alpha*newEnergy) / sqrt(lamb)*(A*log(ab) + B) / C;
    double a = N*c*factor;
    return a;
}

double mainArea::acousticScatteringX(double energy)
{
    double lamb = energy*(1 + nonparabolicityX*energy);
    double j = (sqrt(2)*pow((effectiveMassX*mq), 1.5) * pow(q, 2.5) * kb * T * pow(sigmaX, 2)) / (pow(soundVelocity, 2)*density*pi*pow(h, 4))*sqrt(lamb)*(1 + 2 * nonparabolicityX*energy);
    return j;
}

double mainArea::intervalleyScatteringNVXtoLAbsorbation(double energy)
{
    double N = 1 / (exp(phononXL / (kb*T / q)) - 1);
    double finalMass = effectiveMassL*mq;
    double c = eqValleysL*pow(deformationPotentialConstantLX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononXL)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononXL + energySeparLX;
    double lamb = ef*(1 + nonparabolicityL*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityL*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVXtoLEmission(double energy)
{
    double N = 1 / (exp(phononXL / (kb*T / q)) - 1) + 1;
    double finalMass = effectiveMassL*mq;
    double c = eqValleysL*pow(deformationPotentialConstantLX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononXL)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononXL + energySeparLX;
    double lamb = ef*(1 + nonparabolicityL*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityL*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVXtoGAbsorbation(double energy)
{
    double N = 1 / (exp(phononXG / (kb*T / q)) - 1);
    double finalMass = effectiveMassG*mq;
    double c = eqValleysG*pow(deformationPotentialConstantGX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononXG)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononXG + energySeparGX;
    double lamb = ef*(1 + nonparabolicityG*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityG*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringNVXtoGEmission(double energy)
{
    double N = 1 / (exp(phononXG / (kb*T / q)) - 1)+1;
    double finalMass = effectiveMassG*mq;
    double c = eqValleysG*pow(deformationPotentialConstantGX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononXG)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononXG + energySeparGX;
    double lamb = ef*(1 + nonparabolicityG*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityG*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringEVXtoXEmission(double energy)
{
    if (energy < phononXX)
    {
        return 0.0;
    }
    double N = 1 / (exp(phononXX / (kb*T / q)) - 1) + 1;
    double finalMass = effectiveMassX*mq;
    double c = (eqValleysX - 1)*pow(deformationPotentialConstantXX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononXX)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy - phononXX;
    double lamb = ef*(1 + nonparabolicityX*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityX*ef);
    return ab*factor;
}

double mainArea::intervalleyScatteringEVXtoXAbsorbation(double energy)
{
    double N = 1 / (exp(phononXX / (kb*T / q)) - 1);
    double finalMass = effectiveMassX*mq;
    double c = (eqValleysX - 1)*pow(deformationPotentialConstantXX, 2)*q*sqrt(q) / (sqrt(2)*pi* density *phononXX)*(finalMass / h)*sqrt(finalMass) / h;
    double ab = N*c;
    double ef = energy + phononXX;
    double lamb = ef*(1 + nonparabolicityX*ef);
    double factor = sqrt(lamb)*(1 + 2 * nonparabolicityX*ef);
    return ab*factor;
}

double mainArea::coulombX(double energy)
{
    double alpha = nonparabolicityX;
    double m = mq*effectiveMassX;
    double ld = sqrt(highFrecDielectrConst * kb*T / q * 1 / q * 1 / doppingDensity);
    double ed = h*h / (m * q * 2 * ld*ld);
    double c = doppingDensity*m*q / h*sqrt(2 * m) * 1 / pi * q / h * pow((ld*ld*q / (highFrecDielectrConst*h)), 2) * sqrt(q);
    double a = energy*(1 + alpha*energy);
    double factor = sqrt(a) * (1 + 2 * alpha*energy) * 1 / (1 + 4 * a / ed);
    return c*factor;
}

void mainArea::changeKDirectionCoulombAngle(Electron & e)
{
    double m;
    double alpha;
    switch (e.getValley())
    {
        case 0:
            m = mq*effectiveMassG;
            alpha = nonparabolicityG;
            break;
        case 1:
            m = mq*effectiveMassL;
            alpha = nonparabolicityL;
            break;
        case 2:
            m = mq*effectiveMassX;
            alpha = nonparabolicityX;
            break;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> d(0, 1);

    double enew = e.getEnergy();


    double kxy = sqrt(e.getkX()*e.getkX() + e.getkY()*e.getkY());
    double k = sqrt(kxy*kxy + e.getkZ()*e.getkZ());
    double cth0 = e.getkZ() / k;
    double sth0 = kxy / k;
    double cfi0 = e.getkX() / kxy;
    double sfi0 = e.getkY() / kxy;

    double ge = enew * (1 + alpha * enew);
    double ld = sqrt(highFrecDielectrConst * kb*T / q * 1 / q * 1 / doppingDensity);
    double ed = h*h / (m * q * 2 * ld*ld);
    double rr = d(gen);
    double cth = 1 - 2 * rr / (1 + 4 * (1 - rr)*ge / ed);
    double sth = sqrt(1 - cth*cth);
    double fi = 2 * pi*d(gen);
    double cfi = cos(fi);
    double sfi = sin(fi);
    double kxp = k*sth*cfi;
    double kyp = k*sth*sfi;
    double kzp = k*cth;

    double kx = kxp*cfi0*cth0 - kyp*sfi0 + kzp*cfi0*sth0;
    double ky = kxp*sfi0*cth0 + kyp*cfi0 + kzp*sfi0*sth0;
    double kz = -kxp*sth0 + kzp*cth0;

    e.setKDirection(kx, ky, kz);
}
