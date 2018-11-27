/*
  ScaleSet.cxx
*/
#include <iostream>
#include <string>
#include <cstdio>
#include "AirLattice/ScaleSet.hxx"
#include "AirLattice/Constants.hxx"

ScaleSet* ScaleSet::sInstance=0;

ScaleSet* ScaleSet::get() {
  if (sInstance == 0) {
    sInstance = new ScaleSet();
  }
  return sInstance;
}

ScaleSet::ScaleSet() {
  // Scales
  mL = 1.0;
  mT = 1.0;
  mV = 1.0;
  mRho = 1.0;
  mP = 1.0;
  mTheta0 = 0.0;
  mTheta1 = 0.0;
  mAcceleration = 0.0;
  // Dimensionful parameters
  mKappa = 0.024; // J/(K m s)
  mCp = Constants::Cp_Air; // Cp (air) [J/(K mol)]
  mCp *= 1.0E+3/Constants::MolecularWeight_Air; // [J/(K kg)]
  mRho = 1.1839; // kg/m3
  mNu = 1.5E-5; // m2/s
  mBeta = 0.00369; // K-1
  mG = 9.8; // m/s
  mRe = 1.0;
  mPr = 1.0;
  // Derived parameters
  mParameterA = 1.0;
  mParameterB = 1.0;
}

ScaleSet::~ScaleSet() {
}

double ScaleSet::setScales(double l0, double v0, double rho, 
			   double theta0, double theta1) {
  mL = l0;
  mV = v0;
  mRho = rho;
  mTheta0 = theta0;
  mTheta1 = theta1;
  mP = mRho*mV*mV;
  mT = mL/mV;
  mAcceleration = mV/mT;
  //
  mRe = mL*mV/mNu;
  mPr = mNu*mCp*mRho/mKappa;
  //
  mParameterA = mKappa*mT/(mCp*mRho*mL*mL);
  mParameterB = mBeta*mTheta0*mG/mAcceleration;
}

double ScaleSet::toDimLess(double x, ScaleSet::PropertyType tp) {
}

double ScaleSet::fromDimLess(double x, ScaleSet::PropertyType tp) {
}

double ScaleSet::buoyancyCoefficient() const {
  return mParameterB;
}

double ScaleSet::parameterA() const {
  return mParameterA;
}

void ScaleSet::print() const {
  double kg = mRho*mL*mL*mL;
  double J = mP*mL*mL*mL;

  std::cout << "Scale set" << std::endl;
  std::cout << "# SCALE|PARAMETER          :  VALUE     [UNIT]        Dim.less value" << std::endl;
  std::cout << "SCALES" << std::endl;
  printFormat("L", "m", mL);
  printFormat("V", "m/s", mV);
  printFormat("T", "s", mT);
  printFormat("Acceleration", "m/s^2", mAcceleration);
  printFormat("rho", "kg/m^3", mRho);
  printFormat("P", "J/m^3", mP);
  printFormat("theta0", "K", mTheta0);
  printFormat("theta1", "K", mTheta1);
  std::cout << "Parameters" << std::endl;
  printFormat("g", "m/s^2", mG, mL/(mT*mT) );
  printFormat("kappa", "m^2/s", mKappa, mL*mL/mT);
  printFormat("Cp", "J/K kg", mCp, J/(mTheta1*kg) );
  printFormat("nu", "m^2/s", mNu, mL*mL/mT);
  printFormat("beta0", "1/m^3", mBeta, 1.0/(mL*mL*mL) );
  printFormat("Re", "", mRe, 1.0);
  printFormat("Pr", "", mPr, 1.0);
  std::cout << "Derived parameters" << std::endl;
  printFormat("A (=kappa/(Cp*rho) )", "m^2/s", mParameterA, 1/(mL*mL/mT) );
  printFormat("B (=beta0*theta0*g)", "m/s^2", mBeta*mTheta0*mG, mAcceleration);
}

std::string ScaleSet::format(const std::string& pname, const std::string& unit, 
			     double value_dimful, 
			     double factor_todimful) const {
  char line[200]="";
  char sunit[50]="";

  std::sprintf(sunit, "[%s]", unit.c_str());
  if (factor_todimful>0) {
    double value_dimless = value_dimful/factor_todimful;
    std::sprintf(line, "%-25s: %10.3e %-15s %10.3e", 
		 pname.c_str(), value_dimful, sunit, value_dimless);
  } else {
    std::sprintf(line, "%-25s: %10.3e %-15s", 
		 pname.c_str(), value_dimful, sunit);
  }
  return std::string(line);
}

void ScaleSet::printFormat(const std::string& pname, const std::string& unit, 
			   double value_dimful, double factor_todimful) const {
  std::string s = format(pname, unit, value_dimful, factor_todimful);
  std::cout << "  " << s << std::endl;
}

