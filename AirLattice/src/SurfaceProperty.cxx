/*
  SurfaceProperty.cxx
*/
#include <iostream>
#include <bitset>
#include "AirLattice/SurfaceProperty.hxx"

SurfaceProperty::SurfaceProperty() : mBC() {
}

SurfaceProperty::~SurfaceProperty() {
}

void SurfaceProperty::addBoundaryCondition(BoundaryCondition::Type bctype) {
  mBC.addType(bctype);
}

bool SurfaceProperty::hasBoundaryCondition(BoundaryCondition::Type bctype) const {
  return mBC.hasType(bctype);
}

void SurfaceProperty::print() const {
  std::cout << "SurfaceProperty T=" << T() 
	    << " BoundaryCondition=" << std::bitset<16>(mBC.flags()) << std::endl;
}
