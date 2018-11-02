/*
  LatticeIndexTool.cxx
*/
#include "AirLattice/LatticeIndexTool.hxx"

LatticeIndexTool::LatticeIndexTool() : mEnvironment(0) {
  mNx = 0;
  mNy = 0;
  mNz = 0;
  mNCells = 0;
}

LatticeIndexTool::~LatticeIndexTool() {
}

int LatticeIndexTool::initialize(Environment* env) {
  mEnvironment = env;
  Atmosphere& a = env->atmosphere();
  const int* np = a.NPoints();
  mNx = np[0];
  mNy = np[1];
  mNz = np[2];
  mNCells = mNx*mNy*mNz;
}

int LatticeIndexTool::indexOfCell(int ix, int iy, int iz) const {
  return ( (ix*mNy + iy)*mNz + iz);
}

std::vector<int> LatticeIndexTool::cellIndexAt(int icell) const {
  std::vector<int> v;
  int nyz=0;
  v.assign(3, 0);

  v[0] = icell / (mNy*mNz);
  nyz = icell % (mNy*mNz);
  v[1] = nyz / mNz;
  v[2] = nyz % mNz;
  return v;
}

