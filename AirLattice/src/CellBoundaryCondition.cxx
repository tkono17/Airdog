/*
  CellBoundaryCondition.cxx
*/
#include "AirLattice/CellBoundaryCondition.hxx"


CellBoundaryCondition::CellBoundaryCondition() : mFlags(0) {
}

CellBoundaryCondition::~CellBoundaryCondition() {
}

void CellBoundaryCondition::set(CellBoundaryCondition::Type bctype, bool value) {
  mFlags |= (1<<bctype);
}


bool CellBoundaryCondition::isSet(CellBoundaryCondition::Type bctype) const {
  return (mFlags>>bctype)&0x1;
}

