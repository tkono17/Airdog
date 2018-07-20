/*
  AirProperty.cxx
*/
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/Constants.hxx"

AirProperty::AirProperty() {
  mMolecularWeight = Constants::MolecularWeight_Air;
}

AirProperty::~AirProperty() {
}

void AirProperty::setVelocity(double vx, double vy, double vz) {
  mVelocity[0] = vx;
  mVelocity[1] = vy;
  mVelocity[2] = vz;
}

