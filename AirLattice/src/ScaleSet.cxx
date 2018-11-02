/*
  ScaleSet.cxx
*/
#include "AirLattice/ScaleSet.hxx"

ScaleSet::ScaleSet() {
  mL = 0.0;
  mT = 0.0;
  mV = 0.0;
  mRho = 0.0;
  mP = 0.0;
  mParameterA = 0.0;
  mTheta0 = 0.0;
  mTheta1 = 0.0;
  mAcceleration = 0.0;
}

ScaleSet::~ScaleSet() {
}

double ScaleSet::setScales(double l0, double v0, double rho, 
			   double theta0, double theta1) {
  mL = l0;
  mV = v0;
  mRho = 0.0;
  mTheta0 = theta0;
  mTheta1 = theta1;
  //
  mT = mL/mV;
  mP = mRho*mV*mV;
  mParameterA = mL*mV;
  mAcceleration = mV/mT;
}

double ScaleSet::toDimLess(double x, ScaleSet::PropertyType tp) {
}

double ScaleSet::fromDimLess(double x, ScaleSet::PropertyType tp) {
}

