#include <DUGKSAmr.H>

namespace EoS {
// The chemical potential here is the bulk chemical potential, i.e. \mu_0 = d{\psi}\d{\rho}

//  Van Der Waals Equation Li2021PRE [103, 013304 (2021)]

    const double VdWa = 9./392.;
    const double VdWb = 2./21.;
    const double VdWR = 1.0;
    const double VdWRhoc = 3.5;
    const double VdWPc   = 3./32.;
    const double VdWTc   = 1./14.;

// Redlich-Kwong-Soave Wen2017PRE [95, 063305 (2017)]

    const double RKSa = 2./49.;
    const double RKSb = 2./21.;
    const double RKSR = 1.0;
    const double RKSomega = 0.344;
    const double RKSTc = 0.0868612065393;
    const double RKSPc = 0.0790193768129;

// Peng-Robinson Yuan2006PoF[https://doi.org/10.1063/1.2187070]

    const double PRa = 2./49.;
    const double PRb = 2./21.;
    const double PRc = 1.4142135623731; //sqrt(2)
    const double PRR = 1.0;
    const double PRomega = 0.344;
    const double PRTc = 0.072922;
    const double PRPc = 0.05957;

// Carnahan-Starling  Yuan2006PoF[https://doi.org/10.1063/1.2187070]

    const double CSa = 1.0;
    const double CSb = 4.0;
    const double CSR = 1.0;
    const double CSTc = 0.0943330646786;
    const double CSPc = 0.004416438255591375;

// Pseudo-Potential

    const double PPa = 1.0;
    const double PPG = 4.3;
}

// - Van der Waals
double EoS::PressureVdWIdeal (const double Density, const double Tr)
{
    return Density*Tr*VdWTc*VdWR;
}
double EoS::PressureVdW (const double Density, const double Tr)
{
    return Density*Tr*VdWTc*VdWR/(1.0-VdWb*Density) - VdWa*Density*Density;
}

double EoS::FreeEnergyDensityVdW (const double Density, const double Tr)
{
    return Density*VdWR*Tr*VdWTc*std::log(Density/(1-VdWb*Density)) - VdWa*Density*Density;
}

double EoS::ChemicalPotentialVdW  (const double Density, const double Tr)
{
    return (VdWR*Tr*VdWTc*(std::log(Density/(1.-VdWb*Density)) + 1./(1.-VdWb*Density)) - 2*VdWa*Density);
}

// - Redlich-Kwong-Soave
double EoS::AlphaRKS(const double Tr)
{
    const double tmp = 1. + (0.48+1.574*RKSomega-0.176*RKSomega*RKSomega)*(1.0-std::sqrt(Tr));
    return tmp*tmp;
}

double EoS::PressureRKS (const double Density, const double Tr)
{
    return Density*Tr*RKSTc*RKSR/(1.0-RKSb*Density) - RKSa*AlphaRKS(Tr)*Density*Density/(1.0+RKSb*Density);
}

double EoS::FreeEnergyDensityRKS (const double Density, const double Tr)
{
    return (RKSR*Tr*RKSTc*Density*std::log(Density/(1.-RKSb*Density))
            -RKSa*AlphaRKS(Tr)*Density*std::log(1.+RKSb*Density)/RKSb);
}

double EoS::ChemicalPotentialRKS (const double Density, const double Tr)
{
    const double minusbD = 1.0 - RKSb*Density;
    const double plusbD  = 1.0 + RKSb*Density;

    return (RKSR*Tr*RKSTc*(std::log(Density/minusbD) + 1./minusbD)
    - RKSa*AlphaRKS(Tr)*log(plusbD)/RKSb - RKSa*AlphaRKS(Tr)*Density/plusbD);
}

// - Peng-Robinson

double EoS::AlphaPR(const double Tr)
{
    const double tmp = 1. + (0.37464 + 1.54226*PRomega - 0.26992*PRomega*PRomega)*(1.0-std::sqrt(Tr));
    return tmp*tmp;
}

double EoS::PressurePR (const double Density, const double Tr)
{
    return  Density*Tr*PRTc*PRR/(1.0-PRb*Density)
            - PRa*AlphaPR(Tr)*Density*Density/(1.+(2.-PRb*Density)*PRb*Density);
}

double EoS::FreeEnergyDensityPR (const double Density, const double Tr)
{
    return  PRR*Tr*PRTc*Density*std::log(Density/(1.0-PRb*Density))
            - PRa*AlphaPR(Tr)*Density*0.5*std::log((PRc-1+PRb*Density)/(PRc+1-PRb*Density))/(PRc*PRb);
}
double EoS::ChemicalPotentialPR (const double Density, const double Tr)
{
    return  PRR*Tr*PRTc*std::log(Density/(1.0-PRb*Density)) + PRR*Tr*PRTc/(1.0-PRb*Density)
            - PRa*AlphaPR(Tr)*0.5*std::log((PRc-1+PRb*Density)/(PRc+1-PRb*Density))/(PRc*PRb)
            - PRa*AlphaPR(Tr)*Density/(1.+(2.-PRb*Density)*PRb*Density);
}

// - Carnahan-Starling

double EoS::PressureCS (const double Density, const double Tr)
{
    const double c = 1 - CSb*Density/4.0;
    return Density*CSR*Tr*CSTc*(1+(1-c)*(1+c*(1-c)))/(c*c*c) - CSa*Density*Density;
}

double EoS::FreeEnergyDensityCS (const double Density, const double Tr)
{
    const double c = 1 - CSb*Density/4.0;
    return Density*CSR*Tr*CSTc*((1+2.*c)/(c*c)+std::log(Density)) - CSa*Density*Density;
}

double EoS::ChemicalPotentialCS (const double Density, const double Tr)
{
    const double c = 1 - CSb*Density/4.0;
    return CSR*Tr*CSTc*((2+c)/(c*c*c) + std::log(Density) + 1.) - 2.*CSa*Density;
}

// double well potential
double EoS::ChemicalPotentialDW (const double order, const double orderA, const double orderB)
{
    return (order - orderA)*(order - orderB)*(2*order - orderA - orderB);
}

// - pseudo-potential
double EoS::PseudoPotential(const double Density, const double Tr)
{
    return std::log(-0.5/Density);
}