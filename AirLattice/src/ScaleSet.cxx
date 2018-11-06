/*
  ScaleSet.cxx
*/
#include <iostream>
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
  mL = 1.0;
  mT = 1.0;
  mV = 1.0;
  mRho = 1.0;
  mP = 1.0;
  mTheta0 = 0.0;
  mTheta1 = 0.0;
  mAcceleration = 0.0;
  //
  mKappa = 0.02; // J/(K m s)
  mCp = Constants::Cp_Air; // Cp (air) [J/(K mol)]
  mCp *= 1.0E+3/Constants::MolecularWeight_Air;
  mRho = 1.1839; // kg/m3
  mNu = 0.2E-4; // m2/s
  mBeta = 0.00369; // K-1
  mG = 9.8; // m/s
  mRe = 1.0;
  mPr = 1.0;
}

ScaleSet::~ScaleSet() {
}

double ScaleSet::setScales(double l0, double v0, double rho, 
			   double theta0, double theta1) {
  mL = l0;
  mV = v0;
  mTheta0 = theta0;
  mTheta1 = theta1;
  mP = mRho*mV*mV;
  //
  mT = mL/mV;
  mP = mRho*mV*mV;
  mAcceleration = mV/mT;
  mRe = mL*mV/mNu;
  mPr = mNu*mCp*mRho/mKappa;
}

double ScaleSet::toDimLess(double x, ScaleSet::PropertyType tp) {
}

double ScaleSet::fromDimLess(double x, ScaleSet::PropertyType tp) {
}

double ScaleSet::buoyancyCoefficient() const {
  return mBeta*mTheta1*mG/mAcceleration;
}

void ScaleSet::print() const {
  std::cout << "Scale set" << std::endl;
  std::cout << "  L            : " << mL << std::endl;
  std::cout << "  V            : " << mV << std::endl;
  std::cout << "  T            : " << mT << std::endl;
  std::cout << "  Acceleration : " << mAcceleration << std::endl;
  std::cout << "  g            : " << mG << std::endl;
  std::cout << "  Theta0       : " << mTheta0 << std::endl;
  std::cout << "  Theta1       : " << mTheta1 << std::endl;
  std::cout << "  P            : " << mP << std::endl;
  std::cout << "  kappa        : " << mKappa << std::endl;
  std::cout << "  Cp           : " << mCp << std::endl;
  std::cout << "  rho          : " << mRho << std::endl;
  std::cout << "  nu           : " << mNu << std::endl;
  std::cout << "  beta0        : " << mBeta << std::endl;
  std::cout << "  kappa/(Cp*rho) : " << (mKappa/(mCp*mRho) ) << std::endl;
  std::cout << "  Re           : " << mRe << std::endl;
  std::cout << "  Pr           : " << mPr << std::endl;
  std::cout << "  B            : " << (mBeta*mTheta1*mG)/mAcceleration << std::endl;
}

