/*
  AirProperty.cxx
*/
#include <iostream>
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/Constants.hxx"

AirProperty::AirProperty() {
  mPressure = 1013E+2; // Pa
  mDensity = 1.293E-3; // g/cm3
  mTemperature = 25.0 + Constants::T0; // K
  for (int i=0; i<3; ++i) {
    mVelocity[i] = 0.0; // m/s
  }
  mVaporPressure = 0.0; // Pa
  mWaterDensity = 0.0; // g/cm3

  mMolecularWeight = Constants::MolecularWeight_Air; // g/mol
  mSpecificHeat = Constants::Cp_Air; // Cp (air) [J/(K mol)]
  //  mHeatConductivity = 0.1; // J/(K m2 s)
  mHeatConductivity = 0.02; // J/(K m s)
}

AirProperty::~AirProperty() {
}

void AirProperty::setVelocity(double vx, double vy, double vz) {
  mVelocity[0] = vx;
  mVelocity[1] = vy;
  mVelocity[2] = vz;
}

void AirProperty::setDensityFromPT() {
  double rho_gm3 = mPressure*mMolecularWeight/(Constants::R*mTemperature);
  mDensity = rho_gm3*1.0E-6;
}

double AirProperty::molDensitySI() const {
  double c = 1.0E+6; // 1/cm3 -> 1/m3
  return mDensity*c/mMolecularWeight;
}

double AirProperty::heatCapacity(double volume) const {
  double n = molDensitySI();
  return mSpecificHeat * n * volume;
}

void AirProperty::updateVariable(PropertyType pt, double value) {
  mVariableHolder[pt] = value;
}

void AirProperty::updateVariables() {
  std::map<int, double>::const_iterator p;
  for (p=mVariableHolder.begin(); p!=mVariableHolder.end(); ++p) {
    if (p->first == kT) {
      setTemperature(p->second);
    } else if (p->first == kP) {
      setPressure(p->second);
    } else if (p->first == kRho) {
      setDensity(p->second);
    } else if (p->first == kU) {
      setVelocityU(p->second);
    } else if (p->first == kV) {
      setVelocityV(p->second);
    } else if (p->first == kW) {
      setVelocityW(p->second);
    }
  }
  mVariableHolder.clear();
}

void AirProperty::print() const {
  std::cout << "AirProperty: " 
	    << "P=" << P() << ", T=" << T() << ", rho=" << rho() 
	    << ", (u, v, w)=(" << u() << ", " << v() << ", " << w() << ")" 
	    << std::endl;
}
