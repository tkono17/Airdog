/*
  BoundaryCondition.cxx
*/
#include "AirLattice/BoundaryCondition.hxx"

BoundaryCondition::BoundaryCondition() : mFlags(0) {
}

BoundaryCondition::~BoundaryCondition() {
}

void BoundaryCondition::addType(BoundaryCondition::Type bctype) {
  mFlags |= (1 << bctype);
}

bool 
BoundaryCondition::hasType(BoundaryCondition::Type bctype) const {
  return ( (mFlags >> bctype) & 0x1);
}

